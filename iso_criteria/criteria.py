from fast3tree import fast3tree
import numpy as np



def find_pairs(points, D):
    pairs = []
    with fast3tree(points) as tree:
        for i, p in enumerate(points):
            idx = tree.query_radius(p, D, True, 'index')
            idx = idx[idx > i]
            pairs.extend(((i, j) for j in idx))
    return pairs


def find_periodic_midpoint(p1, p2, box_size):
    ndim = len(p1)
    mid = []
    for i in range(ndim):
        non_pbc_dist = np.abs(p1[i] - p2[i])
        if non_pbc_dist < box_size/2.0:
            mid.append(p1[i] - np.sign(p1[i] - p2[i])*non_pbc_dist/2.0)
        else:
            pbc_dist = box_size - non_pbc_dist
            mid.append((p1[i] + np.sign(p1[i] - p2[i])*pbc_dist/2.0 + box_size)%box_size)
    return mid


def isolated(ps, D, points, masses):
    print 'before: ', len(ps)
    iso_pairs = []
    for pair in ps:
        ptA = points[pair[0]]
        ptB = points[pair[1]]
        mid = find_periodic_midpoint(ptA, ptB, 420)
        with fast3tree(points) as tree:
            idx = tree.query_radius(mid, D, True, 'index')
            masses_in_iso_region = np.array([masses[ind] for ind in idx])
            masses_in_iso_region.sort()
            masses_in_iso_region = masses_in_iso_region[::-1]
            massA = masses[pair[0]]
            massB = masses[pair[1]]
            if ((massA==masses_in_iso_region[0] and massB==masses_in_iso_region[1]) or \
               (massA==masses_in_iso_region[1] and massB==masses_in_iso_region[0])):
                iso_pairs.append(pair)
    print 'after: ', len(iso_pairs)
    return iso_pairs

def label_MW_M31(pairs, masses):
    all_MW = []
    all_M31 = []
    for pair in pairs:
        if masses[pair[0]] > masses[pair[1]]:
            all_M31.append(pair[0])
            all_MW.append(pair[1])
        else:
            all_MW.append(pair[0])
            all_M31.append(pair[1])
    return all_MW, all_M31


def find_M33(all_points, D, all_masses, host_inds, host_points, M31_given):
    all_M33 = []
    if M31_given:
        for M31_ind in host_inds:
            M31_pos = host_points[M31_ind]
            with fast3tree(all_points) as tree:
                idx = tree.query_radius(M31_pos, D, True, 'index')
                M33_masses = [all_masses[x] for x in idx]
                if len(idx)==0:
                    M33 = -1
                else:
                    M33 = idx[M33_masses.index(max(M33_masses))]
                all_M33.append(M33)
    else:
        all_M31 = []
        all_MW = []
        for pair in host_inds:
        M33_pair = []
        M33_mass_pair = []
        for elem in pair:
            M31_pos = host_points[elem]
            with fast3tree(all_points) as tree:
                idx = tree.query_radius(M31_pos, D, True, 'index')
                M33_masses = [all_masses[x] for x in idx]
                if len(idx)==0:
                    M33 = -1
                    mass = -1
                else:
                    mass_ind = M33_masses.index(max(M33_masses))
                    M33 = idx[mass_ind]
                    mass = all_masses[m33_masses[mass_ind]]
                    
                M33_pair.append(M33)
                M33_mass_pair.append(mass)
                if M33_mass_pair[0] > M33_mass_pair[1]:
                    all_M33.append(M33_pair[0])
                    all_M31.append(pair[0])
                    all_MW.append(pair[1])
                else:
                    all_M33.append(M33_pair[1])
                    all_M31.append(pair[1])
                    all_MW.append(pair[0])
    if not M31_given: return all_M33, all_M31, all_MW
    return all_M33
