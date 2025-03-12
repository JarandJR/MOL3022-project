import numpy as np
import pandas as pd
import logomaker
import matplotlib.pyplot as plt
import sys
import json
import datetime

def plot_logo(pwm_json, filename):
    pwm = np.array(json.loads(pwm_json))
    df = pd.DataFrame(pwm.T, columns=["A", "C", "G", "T"])

    crp_logo = logomaker.Logo(df, shade_below=.5, fade_below=.5, font_name='Arial Rounded MT Bold')
    crp_logo.style_spines(visible=False)
    crp_logo.style_spines(spines=['left', 'bottom'], visible=True)
    crp_logo.style_xticks(rotation=90, fmt='%d', anchor=0)
    crp_logo.ax.set_ylabel("Position Weight Matrix Weights", labelpad=-1)
    crp_logo.ax.xaxis.set_ticks_position('none')
    crp_logo.ax.xaxis.set_tick_params(pad=-1)
    
    crp_logo.ax.set_ylim(0, 1)
    
    num_positions = df.shape[0]
    crp_logo.ax.set_xticks(range(num_positions))
    crp_logo.ax.set_xticklabels([str(i+1) for i in range(num_positions)])

    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    filename = f"results/{filename}_sequence_logo_{timestamp}.png"
    
    plt.savefig(filename)
    print(f"Saved the logo as {filename}")
    plt.show()

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python plot_logo.py '<pwm_json>'")
        sys.exit(1)
    pwm_json = sys.argv[1]
    filename = sys.argv[2]
    plot_logo(pwm_json, filename)
