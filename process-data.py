#!/usr/bin/env python3
from sage.all import *
import os, re, collections

colors = {
        'c': 'blue',
        'l': 'green',
    }

################################################################

data = collections.defaultdict(list)

ps, invs = set(), set()

def parse(filename):

    lines = iter(open(f'out/{filename}'))
    variant, _, p, l, inv = re.search(r'\./main big([cl]) (w) ([0-9]+) ([0-9]+) ([0-9]+)', next(lines)).groups()
    p, l, inv = Sequence([p, l, inv], ZZ)

    #TODO check correctness

    ps.add(p)
    invs.add(inv)
    if len(ps) > 1 or len(invs) > 1:
        raise NotImplementedError

    for line in lines:
        if (m := re.search(r'user\s+(.+)', line)):
            break
    else:
        assert False

    time = 0
    t, = m.groups()
    while t:
        num = ''
        while t[0] in '0123456789.':
            num += t[0]
            t = t[1:]
        unit = t[0]
        t = t[1:]
        time += float(num) * {'s':1, 'm':60}[unit]

    return variant, l, time

for _,res in parallel(os.cpu_count())(parse)(os.listdir('out/')):
    variant, l, time = res
    data[variant].append((l, time))

################################################################

lmax = math.ceil(1.15 * max(l for vs in data.values() for l,_ in vs))
tmax = math.ceil(1.35 * max(t for vs in data.values() for _,t in vs))

lstep = 1000
tstep = 86400*10
assert tstep % 86400 == 0
tlabels = [
        f'{t//86400} d'
        for t in range(0, (tmax//tstep+2)*tstep, tstep)
    ]

myplot = scatter_plot([],
                      ticks=[lstep, list(range(0, (tmax//tstep+2)*tstep, tstep))],
                      tick_formatter=[None, tlabels],
                      fontsize=20,
                      gridlines='major',
                      gridlinesstyle={'color': (.9,)*3, 'linestyle': '-'},
              )

myplot.axes_width(2)

def fmttime(t, sep=","):
    if t < 60:
        return fr'{t:7}{sep}s'
    elif t < 3600:
        return fr'{t/60:7.2f}{sep}min'
    elif t < 86400:
        return fr'{t/3600:7.2f}{sep}h'
    else:
        return fr'{t/86400:7.2f}{sep}d'

tottime = 0

for k,vs in sorted(data.items()):

    tottime += sum(t for _,t in vs)

    for l,t in sorted(vs):
        tstr = fmttime(t)
        print(f'{k} {l:9} {tstr}')
    myplot += scatter_plot(vs, figsize=20, marker='.', markersize=400, facecolor=colors[k], edgecolor=None)

    as_ = list(var('a0 a1 a2 a3 a4 a5'))
    x = var('x')

    model = sum(a*x**i for i,a in enumerate(as_))
    fit = find_fit(vs, model, parameters=as_, variables=[x], solution_dict=True)
    print('fit:', [fit[a] for a in as_])
    ff = lambda xx: sum(fit[a]*xx**i for i,a in enumerate(as_))
    myplot += plot(ff, color=colors[k], alpha=.5, linestyle='dotted', thickness=3, xmin=0, xmax=lmax)
    print()

print(f'total time:', tottime / 86400, 'core days')

################################################################

myplot.save('plot.png')

################################################################

# formatting the LaTeX table

print()

ls = sorted({l for k,vs in data.items() for l,_ in vs})
assert set(data.keys()) == {'l','c'}
for l in ls:
    t1 = dict(data['l'])[l]
    t2 = dict(data['c'])[l]
    lstr = f'{l:9}'
    t1str = fr'{fmttime(t1,"\\,"):15}'
    t2str = fr'{fmttime(t2,"\\,"):15}'
    rstr = f'{t1 / t2 : .1f}'
    print(fr'{lstr} & {t1str} & {t2str} & {rstr} \\')

################################################################

