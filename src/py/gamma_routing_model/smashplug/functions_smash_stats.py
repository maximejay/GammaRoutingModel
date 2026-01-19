import numpy as np
import math



def rmse(obs,sim,lacuna=-99.):
    
    mask_lacuna=(obs!=lacuna)
    
    rmse=np.nan
    if isinstance(obs,np.ndarray) and isinstance(sim,np.ndarray):
        
        if (len(obs)==len(sim)):
            rmse=math.sqrt((1./len(mask_lacuna))*np.sum((obs-sim)**2.,where=mask_lacuna))
        else:
            raise ValueError(
                f"Error: len(obs)!=len(sim)"
            )
    else:
        raise ValueError(
                f"Error: obs and sim must be an instance of np.ndarray"
            )
            
    return rmse



def nrmse(obs,sim):
    
    res_rmse=rmse(obs,sim)
    mean_obs=np.mean(obs)
    nrmse=res_rmse/mean_obs
    
    return nrmse


#relative root square error
def rrse(obs,sim,lacuna=-99.):
    
    mask_lacuna=(obs!=lacuna)
    
    rrse=np.nan
    if isinstance(obs,np.ndarray) and isinstance(sim,np.ndarray):
        
        if (len(obs)==len(sim)):
            rrse=math.sqrt(np.sum(((obs-sim)/obs)**2.,where=mask_lacuna))
        else:
            raise ValueError(
                f"Error: len(obs)!=len(sim)"
            )
    else:
        raise ValueError(
                f"Error: obs and sim must be an instance of np.ndarray"
            )
            
    return rrse


def nse(obs,sim,lacuna=-99.):
    
    mask_lacuna=(obs!=lacuna)
    nse=0
    
    if isinstance(obs,np.ndarray) and isinstance(sim,np.ndarray):
        if (len(obs)==len(sim)):
            mean_obs=np.mean(obs,where=mask_lacuna)
            numerator=np.sum((obs-sim)**2.,where=mask_lacuna)
            denominator=np.sum((obs-mean_obs)**2.,where=mask_lacuna)
            if denominator>0.:
                nse=1-numerator/denominator
        else:
            raise ValueError(
                f"Error: len(obs)!=len(sim)"
            )
    else:
        raise ValueError(
                f"Error: obs and sim must be an instance of np.ndarray"
            )
            
    return nse



def nnse(obs,sim):
    
    res_nse=nse(obs,sim)
    nnse=1./(2.-res_nse)
    
    return nnse



def persistance(obs,sim,sim0,echeance=0):
    
    persistance=np.nan
    if isinstance(obs,np.ndarray) and isinstance(sim,np.ndarray):
        if (len(obs)==len(sim)):
            numerator=(sim[echeance]-obs[echeance])**2.
            denominator=(sim0-obs[echeance])**2.
            if denominator>0.:
                persistance=1.-(numerator/denominator)
        else:
                raise ValueError(
                    f"Error: len(obs)!=len(sim)"
                )
    else:
        raise ValueError(
                f"Error: obs and sim must be instance of np.ndarray"
            )
    
    return persistance



def true_persistance(obs,sim,sim0,lacuna=-99.):
    
    mask_lacuna=(obs!=lacuna)
    
    persistance=np.nan
    if isinstance(obs,np.ndarray) and isinstance(sim,np.ndarray):
        if (len(obs)==len(sim)):
            numerator=np.sum((sim[:]-obs[:])**2.,where=mask_lacuna)
            denominator=np.sum((sim0[:]-obs[:])**2.,where=mask_lacuna)
            if denominator>0.:
                persistance=1.-(numerator/denominator)
        else:
                raise ValueError(
                    f"Error: len(obs)!=len(sim)"
                )
    else:
        raise ValueError(
                f"Error: obs and sim must be instance of np.ndarray"
            )
    
    return persistance



def net_persistance(obs,sim,sim0,echeance=0):
    
    persistance=np.nan
    if isinstance(obs,np.ndarray) and isinstance(sim,np.ndarray):
        if (len(obs)==len(sim)):
            numerator=(sim[echeance]-obs[echeance])**2.
            denominator=(sim0-obs[echeance])**2.
            if denominator>0.:
                persistance=(numerator/denominator)
        else:
                raise ValueError(
                    f"Error: len(obs)!=len(sim)"
                )
    else:
        raise ValueError(
                f"Error: obs and sim must be instance of np.ndarray"
            )
    
    return persistance




def gain(obs,sim,sim_ref,echeance=0):
    
    gain=np.nan
    if isinstance(obs,np.ndarray) and isinstance(sim,np.ndarray):
        if (len(obs)==len(sim)):
            numerator=(sim[echeance]-obs[echeance])**2.
            denominator=(sim_ref[echeance]-obs[echeance])**2.
            if denominator>0.:
                gain=1.-(numerator/denominator)
        else:
                raise ValueError(
                    f"Error: len(obs)!=len(sim)"
                )
    else:
        raise ValueError(
                f"Error: obs and sim must be an instance of np.ndarray"
            )
    
    return gain




def true_gain(obs,sim,sim_ref,lacuna=-99.):
    
    mask_lacuna=(obs!=lacuna)
    
    gain=np.nan
    if isinstance(obs,np.ndarray) and isinstance(sim,np.ndarray):
        if (len(obs)==len(sim)):
            numerator=np.sum((sim[:]-obs[:])**2.,where=mask_lacuna)
            denominator=np.sum((sim_ref[:]-obs[:])**2.,where=mask_lacuna)
            if denominator>0.:
                gain=1.-(numerator/denominator)
        else:
                raise ValueError(
                    f"Error: len(obs)!=len(sim)"
                )
    else:
        raise ValueError(
                f"Error: obs and sim must be an instance of np.ndarray"
            )
    
    return gain




def net_gain(obs,sim,sim_ref,echeance=0):
    
    gain=np.nan
    if isinstance(obs,np.ndarray) and isinstance(sim,np.ndarray):
        if (len(obs)==len(sim)):
            numerator=(sim[echeance]-obs[echeance])**2.
            denominator=(sim_ref[echeance]-obs[echeance])**2.
            if denominator>0.:
                gain=(numerator/denominator)
        else:
                raise ValueError(
                    f"Error: len(obs)!=len(sim)"
                )
    else:
        raise ValueError(
                f"Error: obs and sim must be an instance of np.ndarray"
            )
    
    return gain


