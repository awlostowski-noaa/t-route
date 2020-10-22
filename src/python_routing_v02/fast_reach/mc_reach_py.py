import mc_reach as mcr
import numpy as np
import reach
import utils

def __compute_reach_kernel(qup, quc, nreach, input_buf, output_buf):
    """
    Kernel to compute reach.

    Input buffer is array matching following description:
    axis 0 is reach
    axis 1 is inputs in th following order:
        qlat, dt, dx, bw, tw, twcc, n, ncc, cs, s0, qdp, velp, depthp

        qup and quc are initial conditions.

    Output buffer matches the same dimsions as input buffer in axis 0
    Input is nxm (n reaches by m variables)
    Ouput is nx3 (n reaches by 3 return values)
        0: current flow, 1: current depth, 2: current velocity
    """
    
    for i in range(nreach):
        qlat = input_buf[i, 0] # n x 1
        dt = input_buf[i, 1] # n x 1
        dx = input_buf[i, 2] # n x 1
        bw = input_buf[i, 3]
        tw = input_buf[i, 4]
        twcc =input_buf[i, 5]
        n = input_buf[i, 6]
        ncc = input_buf[i, 7]
        cs = input_buf[i, 8]
        s0 = input_buf[i, 9]
        qdp = input_buf[i, 10]
        velp = input_buf[i, 11]
        depthp = input_buf[i, 12]

        out = reach.compute_segment_kernel(
                    dt,
                    qup,
                    quc,
                    qdp,
                    qlat,
                    dx,
                    bw,
                    tw,
                    twcc,
                    n,
                    ncc,
                    cs,
                    s0,
                    velp,
                    depthp,
        )

        output_buf[i, 0] = quc = out["qdc"]
        output_buf[i, 1] = out["velc"]
        output_buf[i, 2] = out["depthc"]

        qup = qdp


def compute_network(nsteps, reaches, connections, 
    parameter_idx, parameter_cols, parameter_values, 
    qlat_values,
    reach_groups,
    reach_group_cache_sizes,
    assume_short_ts=False,
):
    """
    Compute network

    Args:
        nsteps (int): number of time steps
        reaches (list): List of reaches
        connections (dict): Network
        parameter_idx (ndarray): a 1D sorted index for parameter_values
        parameter_values (ndarray): a 2D array of data inputs (nodes x variables)
        qlats (ndarray): a 2D array of qlat values (nodes x nsteps). The index must be shared with parameter_values
        assume_short_ts (bool): Assume short time steps (quc = qup)
        reach_groups (list): number of reaches in each group  
        reach_group_cache_sizes: number of segments in each group

    Notes:
        Array dimensions are checked as a precondition to this method.
    """
    # Check shapes
    if qlat_values.shape[0] != parameter_idx.shape[0] or qlat_values.shape[1] != nsteps:
        raise ValueError(f"Qlat shape is incorrect: expected ({parameter_idx.shape[0], nsteps}), got ({qlat_values.shape[0], qlat_values.shape[1]})")
    if parameter_values.shape[0] != parameter_idx.shape[0] or parameter_values.shape[1] != parameter_cols.shape[0]:
        raise ValueError(f"parameter_values shape mismatch")

    # flowveldepth is 2D float array that holds results
    # columns: flow (qdc), velocity (velc), and depth (depthc) for each timestep
    # rows: indexed by parameter_idx
    flowveldepth = np.zeros((parameter_idx.shape[0], nsteps * 3), dtype='float32')

    # Source columns
    scols = np.array(mcr.column_mapper(parameter_cols), dtype=np.intp)
    
    # hard-coded column. Find a better way to do this
    buf_cols = 13

    # Measure length of all the reaches
    reach_sizes = list(map(len, reaches))
    # For a given reach, get number of upstream nodes
    # cdef list usreach_sizes = [0 for reach in reaches]
    usreach_sizes = [len(connections.get(reach[0], ())) for reach in reaches]

    buf_cache = []

    # reach cache is ordered 1D view of reaches
    # [-len, item, item, item, -len, item, item, -len, item, item, ...]
    reach_cache = np.empty(sum(reach_sizes) + len(reach_sizes), dtype=np.intp)
    # upstream reach cache is ordered 1D view of reaches
    # [-len, item, item, item, -len, item, item, -len, item, item, ...]
    usreach_cache = np.empty(sum(usreach_sizes) + len(usreach_sizes), dtype=np.intp)
    
    ireach_cache = 0
    iusreach_cache = 0
    
    # copy reaches into an array
    for ireach in range(len(reaches)):
        reachlen = reach_sizes[ireach]
        usreachlen = usreach_sizes[ireach]
        reach = reaches[ireach]

        # set the length (must be negative to indicate reach boundary)
        reach_cache[ireach_cache] = -reachlen
        ireach_cache += 1
        bf_results = mcr.binary_find(parameter_idx, reach)
        for bidx in bf_results:
            reach_cache[ireach_cache] = bidx
            ireach_cache += 1

        usreach_cache[iusreach_cache] = -usreachlen
        iusreach_cache += 1
        if usreachlen > 0:
            for bidx in mcr.binary_find(parameter_idx, connections[reach[0]]):
                usreach_cache[iusreach_cache] = bidx
                iusreach_cache += 1
    

    maxreachlen = max(reach_sizes)
    buf = np.empty((maxreachlen, buf_cols), dtype='float32')
    out_buf = np.empty((maxreachlen, 3), dtype='float32')

    drows_tmp = np.arange(maxreachlen, dtype=np.intp)
    timestep = 0

    #with nogil:
    if 1 == 1:
        
        while timestep < nsteps:
            ts_offset = timestep * 3

            ireach_cache = 0
            iusreach_cache = 0
            ireach = 0
            for group_i in range(len(reach_group_cache_sizes)): # loop through prallelizable groupings
                
                ireach_cache_end = ireach_cache + reach_group_cache_sizes[group_i] + reach_groups[group_i]
                
                # continue for all reaches in this groupint - NOTE - we can parallelize this effort.      
                while ireach_cache < ireach_cache_end:
                        
                    reachlen = -reach_cache[ireach_cache]      
                    usreachlen = -usreach_cache[iusreach_cache]
 
                    ireach_cache += 1
                    iusreach_cache += 1
                    #print(ireach_cache, iusreach_cache, np.asarray(reach_cache, dtype=np.intp), np.asarray(usreach_cache, dtype=np.intp))

                    qup = 0.0
                    quc = 0.0     
                    for i in range(usreachlen):
                        
                        quc += flowveldepth[usreach_cache[iusreach_cache + i], ts_offset]
                        if timestep > 0:       
                            qup += flowveldepth[usreach_cache[iusreach_cache + i], ts_offset - 3]

                    buf_view = buf[:reachlen, :]
                    out_view = out_buf[:reachlen, :]
                    drows = drows_tmp[:reachlen]
                    srows = reach_cache[ireach_cache:ireach_cache+reachlen]
                        
                    utils.fill_buffer_column(srows, timestep, drows, 0, qlat_values, buf_view)

                    for i in range(scols.shape[0]):
                            utils.fill_buffer_column(srows, scols[i], drows, i + 1, parameter_values, buf_view)
                    
                    # fill buffer with qdp, depthp, velp
                    if timestep > 0:
                        utils.fill_buffer_column(srows, ts_offset - 3, drows, 10, flowveldepth, buf_view)
                        utils.fill_buffer_column(srows, ts_offset - 2, drows, 11, flowveldepth, buf_view)
                        utils.fill_buffer_column(srows, ts_offset - 1, drows, 12, flowveldepth, buf_view)
                    else:
                        # fill buffer with constant
                        for i in range(drows.shape[0]):
                            buf_view[drows[i], 10] = 0.0
                            buf_view[drows[i], 11] = 0.0
                            buf_view[drows[i], 12] = 0.0

                    if assume_short_ts:
                        quc = qup
                    
                    mcr.compute_reach_kernel(qup, quc, reachlen, buf_view, out_view)
                    
                    # copy out_buf results back to flowdepthvel
                    for i in range(3):
                        utils.fill_buffer_column(drows, i, srows, ts_offset + i, out_view, flowveldepth)
        
                    # Update indexes to point to next reach
                    ireach += 1
                    ireach_cache += reachlen
                    iusreach_cache += usreachlen
                    
            timestep += 1

    return np.asarray(parameter_idx, dtype=np.intp), np.asarray(flowveldepth, dtype='float32')

def compute_network_original(nsteps, reaches, connections, 
    data_idx, data_cols, data_values, 
    qlat_values,
    assume_short_ts=False):
    """
    Compute network
    Args:
        nsteps (int): number of time steps
        reaches (list): List of reaches
        connections (dict): Network
        data_idx (ndarray): a 1D sorted index for data_values
        data_values (ndarray): a 2D array of data inputs (nodes x variables)
        qlats (ndarray): a 2D array of qlat values (nodes x nsteps). The index must be shared with data_values
        assume_short_ts (bool): Assume short time steps (quc = qup)
    Notes:
        Array dimensions are checked as a precondition to this method.
    """
    # Check shapes
    if qlat_values.shape[0] != data_idx.shape[0] or qlat_values.shape[1] != nsteps:
        raise ValueError(f"Qlat shape is incorrect: expected ({data_idx.shape[0], nsteps}), got ({qlat_values.shape[0], qlat_values.shape[1]})")
    if data_values.shape[0] != data_idx.shape[0] or data_values.shape[1] != data_cols.shape[0]:
        raise ValueError(f"data_values shape mismatch")

    # flowveldepth is 2D float array that holds results
    # columns: flow (qdc), velocity (velc), and depth (depthc) for each timestep
    # rows: indexed by data_idx
    flowveldepth = np.zeros((data_idx.shape[0], nsteps * 3), dtype='float32')

    # Source columns
    scols = np.array(mcr.column_mapper(data_cols), dtype=np.intp)
    
    # hard-coded column. Find a better way to do this
    buf_cols = 13

    # Measure length of all the reaches
    reach_sizes = list(map(len, reaches))
    # For a given reach, get number of upstream nodes
    usreach_sizes = [len(connections.get(reach[0], ())) for reach in reaches]

    buf_cache = []

    # reach cache is ordered 1D view of reaches
    # [-len, item, item, item, -len, item, item, -len, item, item, ...]
    reach_cache = np.empty(sum(reach_sizes) + len(reach_sizes), dtype=np.intp)
    # upstream reach cache is ordered 1D view of reaches
    # [-len, item, item, item, -len, item, item, -len, item, item, ...]
    usreach_cache = np.empty(sum(usreach_sizes) + len(usreach_sizes), dtype=np.intp)

    ireach_cache = 0
    iusreach_cache = 0
    # copy reaches into an array
    for ireach in range(len(reaches)):
        reachlen = reach_sizes[ireach]
        usreachlen = usreach_sizes[ireach]
        reach = reaches[ireach]

        # set the length (must be negative to indicate reach boundary)
        reach_cache[ireach_cache] = -reachlen
        ireach_cache += 1
        bf_results = mcr.binary_find(data_idx, reach)
        for bidx in bf_results:
            reach_cache[ireach_cache] = bidx
            ireach_cache += 1

        usreach_cache[iusreach_cache] = -usreachlen
        iusreach_cache += 1
        if usreachlen > 0:
            for bidx in mcr.binary_find(data_idx, connections[reach[0]]):
                usreach_cache[iusreach_cache] = bidx
                iusreach_cache += 1        
    
    maxreachlen = max(reach_sizes)
    buf = np.empty((maxreachlen, buf_cols), dtype='float32')
    out_buf = np.empty((maxreachlen, 3), dtype='float32')

    drows_tmp = np.arange(maxreachlen, dtype=np.intp)
    timestep = 0

    #with nogil:
    if 1 == 1:
        while timestep < nsteps:
            ts_offset = timestep * 3

            ireach_cache = 0
            iusreach_cache = 0
            while ireach_cache < reach_cache.shape[0]:
                
                reachlen = -reach_cache[ireach_cache]
                usreachlen = -usreach_cache[iusreach_cache]

                ireach_cache += 1
                iusreach_cache += 1

                qup = 0.0
                quc = 0.0
                for i in range(usreachlen):
                    quc += flowveldepth[usreach_cache[iusreach_cache + i], ts_offset]
                    if timestep > 0:
                        qup += flowveldepth[usreach_cache[iusreach_cache + i], ts_offset - 3]

                buf_view = buf[:reachlen, :]
                out_view = out_buf[:reachlen, :]
                drows = drows_tmp[:reachlen]
                srows = reach_cache[ireach_cache:ireach_cache+reachlen]

                utils.fill_buffer_column(srows, timestep, drows, 0, qlat_values, buf_view)
                for i in range(scols.shape[0]):
                        utils.fill_buffer_column(srows, scols[i], drows, i + 1, data_values, buf_view)
                    # fill buffer with qdp, depthp, velp
                if timestep > 0:
                    utils.fill_buffer_column(srows, ts_offset - 3, drows, 10, flowveldepth, buf_view)
                    utils.fill_buffer_column(srows, ts_offset - 2, drows, 11, flowveldepth, buf_view)
                    utils.fill_buffer_column(srows, ts_offset - 1, drows, 12, flowveldepth, buf_view)
                else:
                    # fill buffer with constant
                    for i in range(drows.shape[0]):
                        buf_view[drows[i], 10] = 0.0
                        buf_view[drows[i], 11] = 0.0
                        buf_view[drows[i], 12] = 0.0

                if assume_short_ts:
                    quc = qup             
                
                mcr.compute_reach_kernel(qup, quc, reachlen, buf_view, out_view)

                # copy out_buf results back to flowdepthvel
                for i in range(3):
                    utils.fill_buffer_column(drows, i, srows, ts_offset + i, out_view, flowveldepth)

                # Update indexes to point to next reach
                ireach_cache += reachlen
                iusreach_cache += usreachlen
                
            timestep += 1
    return np.asarray(data_idx, dtype=np.intp), np.asarray(flowveldepth, dtype='float32')
