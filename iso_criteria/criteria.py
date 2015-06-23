from fast3tree import fast3tree
import numpy as np


def find_pairs(halos, host_flag, D):
    pairs = []
    points = halos[host_flag][list('xyz')].view(float).reshape(-1, 3)
    with fast3tree(points) as tree:
        for i, p in enumerate(points):
            idx = tree.query_radius(p, D, True, 'index')
            idx = idx[idx > i]
            pairs.extend(((i, j) if halos['mvir'][host_flag[i]] > halos['mvir'][host_flag[j]] else (j, i) for j in idx))
    return pairs

def find_periodic_midpoint(p1, p2, box_size, periodic=True):
    mid = (p1 + p2)*0.5
    if not periodic:
        return mid
    half_box_size = box_size*0.5
    mid[np.fabs(p1 - p2) > half_box_size] += half_box_size
    mid = np.fmod(mid, box_size)
    return mid

def find_largest_sub(tree, center, D, halos, pair_idx, vmax_cut=None, periodic=True):
    idx = tree.query_radius(center, D, periodic, 'index')
    idx = idx[np.in1d(idx, pair_idx, True, True)]
    if len(idx):
        sub_ind = idx[halos['mvir'][idx].argmax()]
        if halos['vmax'][sub_ind] > (vmax_cut or np.inf):
            return -1, 0.0
    else:
        return -1, 0.0
    return sub_ind, halos['mvir'][sub_ind]

def isolated(pairs, halos, host_flag, D_iso, D_M33, box_size, vmax_cut=None, periodic=True):
    print 'before: ', len(pairs)
    larger_host = []
    larger_host_sub = []
    smaller_host = []
    smaller_host_sub = []

    with fast3tree(halos[list('xyz')].view(float).reshape(-1, 3)) as tree:
        for pair in pairs:
            pair_idx = host_flag[list(pair)]
            h1, h2 = halos[pair_idx] #h1 is the larger halo
            p1 = np.fromiter((h1[ax] for ax in 'xyz'), float)
            p2 = np.fromiter((h2[ax] for ax in 'xyz'), float)

            mid = find_periodic_midpoint(p1, p2, box_size, periodic)
            idx = tree.query_radius(mid, D_iso, periodic, 'index')

            biggest = halos['mvir'][idx].argmax()
            idx = np.delete(idx, biggest)
            biggest = halos['mvir'][idx].argmax() #biggest is now second-biggest
            

            # the second index has a smaller mass
            if h2['mvir'] == halos['mvir'][idx[biggest]]:
                larger_host.append(pair[0])
                smaller_host.append(pair[1])
            else:
                continue #no need to search M33

            #find M33
            M33_ind, M33_mass = find_largest_sub(tree, p1, D_M33, halos, pair_idx, vmax_cut, periodic)
            LMC_ind, LMC_mass = find_largest_sub(tree, p2, D_M33, halos, pair_idx, vmax_cut, periodic)
            larger_host_sub.append(M33_ind)
            smaller_host_sub.append(LMC_ind)

    print 'after: ', len(larger_host)
    return larger_host, larger_host_sub, smaller_host, smaller_host_sub


def get_data(halos, host_flag, MW_inds, M31_inds, M33_inds, LMC_inds):
    fn = ['MW_'+nm for nm in halos.dtype.names]
    fn1 = ['M31_'+nm for nm in halos.dtype.names]
    fn2 = ['M33_'+nm for nm in halos.dtype.names]
    fn3 = ['LMC_' +nm for nm in halos.dtype.names]
    fields = np.hstack((np.array(fn), np.array(fn1), np.array(fn2), np.array(fn3)))
    dtype = np.dtype([(f, int if f.endswith('id') else float) for f in fields])
    num_sys = len(MW_inds)
    data = np.empty(num_sys, dtype)
    no_sub_data = np.copy(halos[0])
    for f in no_sub_data.dtype.names: no_sub_data[f]=-1
    mw_dat = halos[host_flag[MW_inds]]
    m31_dat = halos[host_flag[M31_inds]]
    m33_dat = np.array([halos[i] if i>0 else no_sub_data for i in M33_inds])
    lmc_dat = np.array([halos[i] if i>0 else no_sub_data for i in LMC_inds])
    for i in range(len(data)):
        data[i] = tuple(mw_dat[i])+tuple(m31_dat[i])+tuple(m33_dat[i])+tuple(lmc_dat[i])
    return data

def get_trip_data(data):
    trip_data = data[data['M33_id']!=-1]
    return trip_data

def get_quad_data(data):
    quad_data = data[data['M33_id']!=-1]
    quad_data = quad_data[quad_data['LMC_id']!=-1]
    return quad_data

