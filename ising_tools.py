import contextlib
import random

@contextlib.contextmanager
def train_test_configs(bins, n_bins, train_frac):
    """Determine whether the current configuration will be used for 
    training or testing.
    """
    if bins < (train_frac * n_bins):
        spin_confgs_file = open('Xtrain.txt', 'a')
        phases_file = open('ytrain.txt', 'a')
    else:
        spin_confgs_file = open('Xtest.txt', 'a')
        phases_file = open('ytest.txt', 'a')
        
    yield [spin_confgs_file, phases_file]
    
    spin_confgs_file.close()
    phases_file.close()


def write_configs(bins, n_bins, train_frac, temperature, spin_model):
    """Write configurations in a text file.
    
    This model has only two phases 0 or 1 depending if 
    we are above or below the critical temperature of `spin_model`.
        
    Args:
        bins        (int):    Current measurement.
        n_bins      (int):    Total number of measurement bins.          
        train_frac  (float):  Fraction of data to be used for training.
        temperature (float):  Value of temperature.
        spin_model  (class):  An Ising2D instance.
    """
    with train_test_configs(bins, n_bins, train_frac) as files:
        # Multiply the current configuration by either +1 or -1 to ensure 
        # we generate configurations with both positive and negative magnetization
        flip = 2*random.randint(0, 1) - 1
        
        # Loop to write each spin to a single line of the spin_confgs_file
        for index in range(spin_model.n_spins):
            
            current_spin = flip * spin_model.spins[index]    
            # Replace -1 with 0 (to be consistent with the desired format)
            if current_spin == -1:
                current_spin = 0
                
            files[0].write(f'{current_spin}  ')
        
        files[0].write('\n')
        
        phase = 0
        if temperature > spin_model.crit_temperature:
            phase = 1
        files[1].write(f'{phase} \n')