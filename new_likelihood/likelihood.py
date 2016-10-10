from scipy.stats import norm, rice


def _vector_norm(x):
    return np.sqrt((x*x).sum(axis=0))


def _positive_normal(loc=0.0, scale=1.0, size=1):
    if loc < scale:
        raise ValueError("Are you sure you want to do this? Probably not...")
    size = int(size)
    x = np.random.normal(loc, scale, size=size)
    flag = (x <= 0.0)
    k = np.count_nonzero(flag)
    while k:
        x[flag] = np.random.normal(loc, scale, size=k)
        flag = (x <= 0.0)
        k = np.count_nonzero(flag)
    return x


def _kpc2dm(d):
    mu = np.log10(d)
    mu += 2.0
    mu *= 5.0
    return mu


_yr2sec = 31556926.0
_kpc2km = 3.08567758e16
_rad2mas = 206265.0e3
_masyr2kmskpc = _kpc2km/_rad2mas/_yr2sec


def sample_solar_system(R0=(8.29, 0.16), V_R0=(30.2, 0.2), U=(11.1, 1.23), W=(7.25, 0.62), size=1):
    """
    R0 in kpc
    V_RO in km/s/kpc
    U in km/s
    W in km/s
    """
    size = int(size)
    
    if size < 1:
        # return mean value in GC, no random
        pos = np.array([R0[0], 0.0, 0.0], dtype=float)
        vel = np.array([-U[0], V_R0[0]*R0[0], W[0]], dtype=float)
        return pos, vel
        
    pos = np.random.randn(3, size)
    pos /= _vector_norm(pos)
    R = _positive_normal(*R0, size=size)
    
    vel = pos * np.random.normal(-U[0], U[1], size=size)
    
    v = np.cross(np.random.randn(3, size), pos, 0, 0, 0)
    v /= _vector_norm(v)
    w = np.cross(v, pos, 0, 0, 0)
    w /= _vector_norm(w)
    
    v *= np.random.normal(*V_R0, size=size)
    v *= R
    vel += v
    w *= np.random.normal(*W, size=size)
    vel += w
    
    pos *= R
    
    return pos.T, vel.T


def convert_gc_to_hc(gc_pos, gc_vel, solar_pos, solar_vel):
    """
    pos in kpc
    vel in km/s
    """
    pos = gc_pos - solar_pos
    vel = gc_vel - solar_vel
    D = np.sqrt((pos*pos).sum(axis=-1))
    dm_theory = _kpc2dm(D)
    vlos_theory = (vel*pos).sum(axis=-1)
    vlos_theory /= D
    mu_theory = np.sqrt((vel*vel).sum(axis=-1) - vlos_theory*vlos_theory)
    mu_theory /= D
    return dm_theory, vlos_theory, mu_theory
    

def _log_sum(x):
    x = np.asarray(x)
    x = x[np.isfinite(x)]
    if not len(x):
        return -np.inf
    i = x.argmax()
    x_max = x[i]
    x -= x_max
    x = np.exp(x, out=x)
    x[i] = 0.0
    return x_max + np.log1p(x.sum())


def calc_log_likelihood_kinematics(theory, obs):
    score = norm.logpdf(obs['dm'][0], loc=theory['dm'], scale=obs['dm'][1])
    score += norm.logpdf(obs['vlos'][0], loc=theory['vlos'], scale=obs['vlos'][1])
    score += rice.logpdf(obs['mu'][0], b=theory['mu']/obs['mu'][1], scale=obs['mu'][1])
    return _log_sum(score) - np.log(float(len(score)))


_obs = {'LMC':{'vmax': (91.7, 18.8),
              'dm': (18.5, 0.1), #Freedman 2001
              'vlos': (262.2, 3.4), #van der Marel 2002
              'mu': tuple(np.array([np.sqrt((1.910*np.cos(69.19/180.0*np.pi))**2.0+0.229**2.0), 0.047]) * _masyr2kmskpc)}, #Kallivayalil 2013
       'SMC':{'vmax': (65.0, 15.0),
              'dm': (18.99, 0.1), #Cioni 2000
              'vlos': (145.6, 0.6), #Harris 2006
              'mu': tuple(np.array([np.sqrt((0.772*np.cos(72.42/180.0*np.pi))**2.0+1.117**2.0), 0.063]) * _masyr2kmskpc)}, #Kallivayalil 2013
       }


def calc_log_likelihood(host, obs1=_obs['LMC'], obs2=_obs['SMC'], solar_sample_size=100000):
    solar = sample_solar_system(size=solar_sample_size)
    s1hc = dict(zip(('dm', 'vlos', 'mu'), convert_gc_to_hc(host['s1_pos'], host['s1_vel'], *solar)))
    s2hc = dict(zip(('dm', 'vlos', 'mu'), convert_gc_to_hc(host['s2_pos'], host['s2_vel'], *solar)))
    del solar
    
    score1 = calc_log_likelihood_kinematics(s1hc, obs1)
    score1 += norm.logpdf(obs1['vmax'][0], loc=host['s1_vmax'], scale=obs1['vmax'][1])
    
    score1 += calc_log_likelihood_kinematics(s2hc, obs2)
    score1 += norm.logpdf(obs2['vmax'][0], loc=host['s2_vmax'], scale=obs2['vmax'][1])
    
    score2 = calc_log_likelihood_kinematics(s1hc, obs2)
    score2 += norm.logpdf(obs2['vmax'][0], loc=host['s1_vmax'], scale=obs2['vmax'][1])
    
    score2 += calc_log_likelihood_kinematics(s2hc, obs1)
    score2 += norm.logpdf(obs1['vmax'][0], loc=host['s2_vmax'], scale=obs1['vmax'][1])
    
    return _log_sum([score1, score2])


if __name__ == "__main__":
    print convert_gc_to_hc(np.array([-0.378, 0.613, -0.283])*1.0e3, np.array([66.213, -76.264, 45.008]), *sample_solar_system(size=-1))

    print _kpc2dm(770.0), -310.0, np.sqrt((-125.2)**2.0 + (-73.8)**2.0) / 770.0

    print convert_gc_to_hc(np.array([-0.001, -0.041, 0.028])*1.0e3, np.array([-56.519, -224.889, 220.139]), *sample_solar_system(size=-1))

    print _obs['LMC']['dm'][0], _obs['LMC']['vlos'][0], _obs['LMC']['mu'][0]
