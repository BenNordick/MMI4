import _loadorder as _
import glob
from io import StringIO
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import pathlib
import roadrunner
import rrplugins
import warnings

auto = rrplugins.Plugin('tel_auto2000')
last_bif_diagram = None
last_runner = None

def run_bifurcation(runner: roadrunner.RoadRunner, param_name, param_limits, **kwargs):
    global last_bif_diagram, last_runner
    need_fort8 = len(runner.fs()) > 6
    if need_fort8:
        warnings.warn('Extra columns NYI')
    auto.setProperty('SBML', runner.getCurrentSBML())
    auto.setProperty('PreSimulation', False)
    auto.setProperty('ScanDirection', 'Positive')
    auto.setProperty('PrincipalContinuationParameter', param_name)
    auto.setProperty('RL0', param_limits[0])
    auto.setProperty('RL1', param_limits[1])
    auto.setProperty('NMX', 10000)
    auto.setProperty('NPR', 2 if need_fort8 else 50)
    auto.setProperty('DSMIN', 1e-5)
    auto.setProperty('DSMAX', 0.01)
    auto.setProperty('ITNW', 10)
    auto.setProperty('KeepTempFiles', need_fort8)
    for prop, value in kwargs.items():
        auto.setProperty(prop, value)
    tel_fail = False
    tmpdir = pathlib.Path(auto.getProperty('TempFolder'))
    try:
        auto.execute()
        fort7_data = auto.BifurcationDiagram
    except Exception as e:
        if 'StringList' in str(e) or 'Tellurium exception' in str(e):
            tel_fail = True
            fort7_path = str(tmpdir / 'fort.7')
            with open(fort7_path) as fort7:
                lines = []
                for line in fort7:
                    if len(line) > 14:
                        lines.append(f'{line[:14]} {line[14:]}')
                    else:
                        lines.append(line)
                fort7_data = ''.join(lines)
        else:
            raise
    prolog_rows = 11
    if '\n   2' in fort7_data:
        tables = fort7_data.split('\n   0\n')
        fort7_data = tables[1]
        prolog_rows = 0
        warnings.warn(f'Ignoring {len(tables) - 2} extra branch(es)')
    bif_diagram = pd.read_csv(StringIO(fort7_data), delim_whitespace=True, skiprows=prolog_rows, skipinitialspace=True)
    bif_diagram['stable'] = -np.sign(bif_diagram['PT'])
    bif_diagram['row'] = np.abs(bif_diagram['PT'])
    if need_fort8:
        points_data = {}
        current_point_fs = None
        with open(str(tmpdir / 'fort.8')) as fort8:
            for line in fort8:
                line = line.rstrip()
                if line == '': # End of file
                    break
                elif line.startswith('      '): # Point data
                    for part in line.strip().split():
                        if len(current_point_fs) <= len(runner.fs()):
                            current_point_fs.append(float(part))
                else: # Point header
                    parts = line.strip().split()
                    point_id = int(parts[3])
                    current_point_fs = []
                    points_data[point_id] = current_point_fs
        bif_diagram = bif_diagram[bif_diagram['LAB'] > 0].reset_index(drop=True)
        labels = list(bif_diagram['LAB'])
        for u in range(7, len(current_point_fs)):
            u_data = [points_data[i][u] for i in labels]
            bif_diagram[f'U({u})'] = u_data
        if not ('KeepTempFiles' in kwargs and kwargs['KeepTempFiles']):
            for filename in glob.iglob(str(tmpdir / 'fort.*')):
                os.remove(filename)
    last_bif_diagram = bif_diagram
    last_runner = runner
    return bif_diagram

def segment_ranges(bif_diagram):
    segment_last_rows = np.argwhere(np.diff(bif_diagram['stable'])).flatten()
    for start, end in zip([0, *segment_last_rows], [*segment_last_rows, len(bif_diagram) - 1]):
        if end - start > 0:
            yield start, end

def plot_bifurcation_diagram(runner, bif_diagram, ax=None, species=None):
    fig = None
    if ax is None:
        fig, ax = plt.subplots()
    actual_order = runner if isinstance(runner, list) else runner.fs()
    needs_label = True
    for start, end in segment_ranges(bif_diagram):
        segment = bif_diagram.iloc[start:(end + 1)]
        species_ids = list(enumerate(actual_order))
        stable = bif_diagram['stable'][end] > 0
        for i, fs in sorted(species_ids, key=lambda nf: nf[1]):
            if species is not None and fs not in species:
                continue
            u = f'U({i + 1})'
            if u not in bif_diagram.columns:
                continue
            legend = fs if stable and needs_label else None
            if stable:
                color = species[fs] if species and fs in species else f'C{i}'
                ax.plot(segment['PAR(0)'], segment[u], color=color, label=legend)
            else:
                ax.plot(segment['PAR(0)'], segment[u], '--', color='gray', label=legend)
        if stable:
            needs_label = False
    ax.legend()
    return fig, ax

def plot_last_bifurcation(*args, actual_order=None, **kwargs):
    runner_or_order = last_runner if actual_order is None else actual_order
    fig, ax = plot_bifurcation_diagram(runner_or_order, last_bif_diagram, *args, **kwargs)
    ax.set_xlabel(auto.getProperty('PrincipalContinuationParameter'))
    ax.set_ylabel('Concentration (AU)')
    ax.set_xlim(left=auto.getProperty('RL0'), right=auto.getProperty('RL1'))
    ax.set_ylim(bottom=0)
    return fig, ax

def plot_split_bifurcation_diagram(runner_or_order, bif_diagram, species, axs, state_colors):
    actual_order = runner_or_order if isinstance(runner_or_order, list) else runner_or_order.fs()
    next_color_id = 0
    for start, end in segment_ranges(bif_diagram):
        segment = bif_diagram.iloc[start:(end + 1)]
        species_ids = list(enumerate(actual_order))
        stable = bif_diagram['stable'][end] > 0
        if stable:
            color = state_colors[next_color_id]
            next_color_id += 1
        for i, fs in sorted(species_ids, key=lambda nf: nf[1]):
            if species is not None and fs not in species:
                continue
            u = f'U({i + 1})'
            if u not in bif_diagram.columns:
                continue
            ax = axs[species.index(fs)]
            if stable:
                ax.plot(segment['PAR(0)'], segment[u], color=color)
            else:
                ax.plot(segment['PAR(0)'], segment[u], '--', color='gray')

def region_attractor_counts(bif_diagram):
    stable_segments = [(start, end) for start, end in segment_ranges(bif_diagram) if bif_diagram['stable'][end] > 0]
    stops = [(bif_diagram['PAR(0)'][start], 1) for start, _ in stable_segments] + [(bif_diagram['PAR(0)'][end], -1) for _, end in stable_segments]
    results = []
    prev_stop = None
    current_attractors = 0
    for stop, diff in sorted(stops, key=lambda sd: sd[0]):
        if stop != prev_stop and prev_stop is not None:
            results.append((prev_stop, stop, current_attractors))
        prev_stop = stop
        current_attractors += diff
    if all(attractors <= 0 for _, _, attractors in results):
        return [(start, end, abs(attractors)) for start, end, attractors in results]
    else:
        return results
