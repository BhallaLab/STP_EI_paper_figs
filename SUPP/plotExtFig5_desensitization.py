import numpy as np
import matplotlib.pyplot as plt
import argparse

tauCell = 0.01 # 10 ms
ChR2chanOpenTime = 0.001    # Note that this balances out tauChR2chan
tauChR2chan = 0.005 # 5 ms to charge cell when ChR2 chan is open.
tauChR2recovery = 1.5   # 1.5 sec for ChR2 to recover from desensitization
ChR2decrement = 0.0004
Erest = 0
EChR2 = 60
simDt = 0.0005

def calcVm( events, dt, ChR2_basal_desensitization ):
    ret = []
    lastStim = -10.0    # Long time since last event.
    dy = ChR2_basal_desensitization
    Vm = Erest
    G = 1.0
    for idx, rr in enumerate( events ):
        t = idx * dt
        if rr:
            dy =  ChR2_basal_desensitization + ChR2decrement/(t-lastStim)
            #G -= dt * G*( 1-np.exp(-dy) )
            G -= G*( 1-np.exp(-dy) )
            Vm += ChR2chanOpenTime * ( EChR2-Vm ) * G / tauChR2chan
            lastStim = t
        G += dt * ( 1-G )/tauChR2recovery
        Vm += dt * (Erest-Vm)/tauCell
        ret.append( Vm )
    return np.array( ret )

def simulate_A(ax, idx, tau, scale, freq, duration=5.0, dt=simDt):
    """
    Simulate and plot the evolution of A over time for set freq interval.
    
    Parameters:
    tau (float): Time constant.
    stimMax (float): Maximum stimulus amplitude.
    freq (list of float): List of intervals between uniform stimulus events (s).
    duration (float): Total duration of the simulation (s). Default is 10.0.
    dt (float): Time step (s). Default is 0.001.
    """
    
    time = np.arange(0, duration, dt)
    stim_events = np.zeros_like(time)
    event_times = np.arange(0, 32) / freq + 1.0
    event_times = np.insert( event_times, 0, 0.2 )
    event_times = np.insert( event_times, 0, 0.205 )
    event_times = np.insert( event_times, 0, 0.215 )
    event_times = np.insert( event_times, 0, 0.230 )
    event_times = np.insert( event_times, 0, 0.250 )
    event_times = np.insert( event_times, 0, 0.275 )
    event_times = np.insert( event_times, 0, 0.305 )
    event_times = np.insert( event_times, 0, 0.700 )
    stim_indices = (event_times / dt).astype(int)
    stim_events[stim_indices] = 1
    BasalDesensitization = 0.01
    #A = desensitization( stim_events, dt, BasalDesensitization, scale, tau )
    Vm = calcVm( stim_events, dt, 0.01 )

    ax.plot(time, Vm, label='$Vm(t)$', color='blue')
    #ax.plot(time, stim_events * 0.2 + 1.2, label='Stim', color='red')
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('$Vm$ (mV)')
    ax.set_xlim( 0, 1.0 + 40/freq + 2*tauCell )
    #ax.set_ylim( 0, 1.5 )
    ax.set_title(f'Freq = {freq} Hz')
    #ax.legend()


def main():
    parser = argparse.ArgumentParser(description='Simulate A over time with given parameters.')
    parser.add_argument('-t', '--tau', type=float, help='Time constant.', default = 0.4)
    parser.add_argument('-s', '--scale', type=float, help='Scaling for decrement.', default = 0.004)

    args = parser.parse_args()
    #event_intervals = [0.125, 0.05, 0.02]  # List of event intervals to simulate
    freqs = [8, 20, 50]  # List of frequencies
    plt.rcParams.update( {"font.size": 16} )
    fig, axs = plt.subplots(len(freqs), 1, figsize=(8, 4 * len(freqs)))
    for idx, ff in enumerate(freqs):

        simulate_A(axs[idx], idx, tau=args.tau, scale=args.scale, freq=ff)
    plt.tight_layout()
    plt.show()

if __name__ == '__main__':
    main()

