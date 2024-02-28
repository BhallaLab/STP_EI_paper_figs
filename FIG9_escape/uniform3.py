import numpy as np
import matplotlib.pyplot as plt
import argparse

def simulate_A(ax, idx, tau, scale, event_interval, duration=5.0, dt=0.001):
    """
    Simulate and plot the evolution of A over time for set event interval.
    
    Parameters:
    tau (float): Time constant.
    stimMax (float): Maximum stimulus amplitude.
    event_intervals (list of float): List of intervals between uniform stimulus events (s).
    duration (float): Total duration of the simulation (s). Default is 10.0.
    dt (float): Time step (s). Default is 0.001.
    """
    
    time = np.arange(0, duration, dt)
    A = np.ones_like(time)
    stim_events = np.zeros_like(time)
    event_times = np.arange(0, 32) * event_interval + 0.5
    event_times = np.insert( event_times, 0, 0.2 )
    stim_indices = (event_times / dt).astype(int)
    stim_events[stim_indices] = 1
    BasalDesensitization = 0.01
    dy = BasalDesensitization
    lastStim = -10.0 # long time since last event.

    for i in range(1, len(time)):
        stim = stim_events[i-1]
        if stim:
            interval = time[i-1] - lastStim
            dy = BasalDesensitization + scale/interval
            lastStim = time[i-1]
            print( i, lastStim, interval )
        A[i] = A[i-1] + (1 - A[i-1]) * dt / tau - stim*A[i-1] * (1 - np.exp(-dy))

    ax.plot(time, A, label='$A(t)$', color='blue')
    ax.plot(time, stim_events * 0.2 + 1.2, label='Stim', color='red')
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('$A$')
    ax.set_xlim( 0, 0.5 + event_interval * 33 )
    ax.set_ylim( 0, 1.5 )
    ax.set_title(f'Event Interval = {event_interval}s')
    ax.legend()


def main():
    parser = argparse.ArgumentParser(description='Simulate A over time with given parameters.')
    parser.add_argument('-t', '--tau', type=float, help='Time constant.', default = 0.4)
    parser.add_argument('-s', '--scale', type=float, help='Scaling for decrement.', default = 0.004)

    args = parser.parse_args()
    event_intervals = [0.125, 0.05, 0.02]  # List of event intervals to simulate
    plt.rcParams.update( {"font.size": 16} )
    fig, axs = plt.subplots(len(event_intervals), 1, figsize=(8, 4 * len(event_intervals)))
    for idx, event_interval in enumerate(event_intervals):

        simulate_A(axs[idx], idx, tau=args.tau, scale=args.scale, event_interval=event_interval)
    plt.tight_layout()
    plt.show()

if __name__ == '__main__':
    main()

