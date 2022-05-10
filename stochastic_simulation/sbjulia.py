import argparse
import collections
import re

parser = argparse.ArgumentParser()
parser.add_argument('input', type=str, help='input Antimony or Python file')
parser.add_argument('output', type=str, help='output Julia file')
parser.add_argument('--noisespecies', type=str, default='._', help='regex for noisy species')
parser.add_argument('--nonnegativespecies', type=str, default='', help='regex for nonnegative species')
args = parser.parse_args()

is_python = args.input.endswith('.py')
in_sb = not is_python
lines = []
with open(args.input) as f:
    for line in f:
        line = line.rstrip()
        if in_sb:
            if is_python and line.startswith("'''"):
                break
            lines.append(line)
        elif line.endswith(".loada('''"):
            in_sb = True

def is_comment(line):
    return re.match('\\s*//', line) is not None

lines = [line for line in lines if (line.strip() != '' and not is_comment(line))]

def is_reaction(line):
    return re.match('\\w+:', line) is not None

species_reactions = collections.defaultdict(list)
def add_rates(reaction_side, reaction, stoich_multiplier):
    for species_spec in reaction_side.split('+'):
        m = re.match('((\\d+)\\s)?(\\w+)', species_spec.strip())
        if m is not None:
            _, stoich, species = m.groups()
            if stoich is None:
                stoich = 1
            species_reactions[species].append((stoich * stoich_multiplier, reaction))

reaction_lines = [line for line in lines if is_reaction(line)]
reactions = {}
for line in reaction_lines:
    m = re.match('(\\w+):\\s*(.*)\\s+->\\s*(.*);\\s*(.*)', line)
    reaction, reactants, products, ratelaw = m.groups()
    reactions[reaction] = ratelaw
    add_rates(reactants, reaction, -1)
    add_rates(products, reaction, 1)

nonreaction_lines = [line for line in lines if not is_reaction(line)]
param_lines = []
raw_blocks = []
in_raw_block = False
for line in nonreaction_lines:
    if in_raw_block:
        raw_blocks[-1] += '\n' + line
        if line.strip() == 'end':
            in_raw_block = False
    elif line.startswith('function'):
        in_raw_block = True
        raw_blocks.append(line)
    else:
        param_lines.append(line)

definition_replacements = {}
free_parameters = []
parameter_defaults = {}
for line in param_lines:
    for assignment in line.split(';'):
        m = re.match('(\\w*) := (.*)', assignment.strip())
        if m is not None:
            defined_param, expression = m.groups()
            definition_replacements[defined_param] = expression
        else:
            m = re.match('(\\w*)\\s*=\\s*(.*)', assignment.strip())
            if m is not None:
                param, default = m.groups()
                if param not in species_reactions:
                    free_parameters.append(m.groups()[0])
                    parameter_defaults[param] = default

noisy_species = [s for s in species_reactions.keys() if re.match(args.noisespecies, s) is not None]
all_parameters = list(free_parameters)
all_parameters.extend([f'σ{s}' for s in noisy_species])
all_parameters.append('τ')
all_species = list(species_reactions.keys())
species_noises = {}
noise_original_species = {}
for species in noisy_species:
    noise_process = f'ξ{species}'
    species_noises[species] = noise_process
    noise_original_species[noise_process] = species
    all_species.append(noise_process)

with open(args.output, 'w') as f:
    for block in raw_blocks:
        f.write(block)
        f.write('\n\n')
    f.write('function model(du, u, p, t)\n')
    param_unpack = ', '.join(all_parameters)
    f.write(f'    {param_unpack} = p\n')
    species_unpack = ', '.join(all_species)
    f.write(f'    {species_unpack} = u\n')
    nonnegative_species = [s for s in species_reactions.keys() if re.match(args.nonnegativespecies, s) is not None]
    for species in nonnegative_species:
        f.write(f'    {species} = max({species}, 0.)\n')
    f.write('\n')
    for (definition, value) in definition_replacements.items():
        f.write(f'    {definition} = {value}\n')
    f.write('\n')
    for (reaction, ratelaw) in reactions.items():
        f.write(f'    {reaction} = {ratelaw}\n')
    f.write('\n')
    for (s, species) in enumerate(all_species):
        f.write(f'    du[{s + 1}] = ')
        if species in noise_original_species:
            f.write(f'-u[{s + 1}] / τ\n')
        else:
            terms = [f'({stoich})*{reaction}' for (stoich, reaction) in species_reactions[species]]
            if species in species_noises:
                terms.append(species_noises[species])
            total_ratelaw = ' + '.join(terms)
            f.write(f'{total_ratelaw}\n')
    f.write('end\n\n')
    f.write('function stoch(du, u, p, t)\n')
    f.write(f'    {species_unpack} = u\n') # TEMP?
    f.write(f'    {param_unpack} = p\n\n')
    for (s, species) in enumerate(all_species):
        f.write(f'    du[{s + 1}] = ')
        if species in noise_original_species:
            noise_param = 'σ' + noise_original_species[species]
            f.write(f'{noise_param} * sqrt(2.0 / τ) * {noise_original_species[species]}\n') # TEMP?
        else:
            f.write('0\n')
    f.write('end\n\n')
    f.write('species_ids = Dict(\n')
    for (s, species) in enumerate(all_species):
        f.write(f'    "{species}" => {s + 1},\n')
    f.write(')\n\n')
    f.write('param_ids = Dict(\n')
    for (p, param) in enumerate(all_parameters):
        f.write(f'    "{param}" => {p + 1},\n')
    f.write(')\n\n')
    f.write('p_default = [\n')
    for (p, param) in enumerate(all_parameters):
        if param in parameter_defaults:
            f.write(f'    {parameter_defaults[param]}')
        elif param == 'τ':
            f.write('    0.1')
        else:
            f.write('    0')
        f.write(f', # {p + 1} {param}\n')
    f.write(']\n')
