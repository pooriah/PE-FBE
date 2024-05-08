# Investigating Bubble-Induced Overpotential, Current Non-Uniformity, and Gas Cross-over in Flow-based Water Electrolyzers: A Numerical Study

## Build

**Aphros** can be installed from the following link:

[**Aphros** code: https://github.com/cselab/aphros](https://github.com/cselab/aphros)

Before compiling the code, follow the instructions below.

1. The following part of **Aphros** has been modified:

path: `src/solver/electro.ipp`

changes in `Electro<EB_>::Imp` -> `void Step(...`

this line:

```
t.ff_resist[cf] = 1 / (1 / r2 * ff_vf[cf] + 1 / r1 * (1 - ff_vf[cf]));
```

has been changed to:

```
	if(ff_vf[cf]>1)
        {
                Scal temp=1-1e-6;
                t.ff_resist[cf] = 1 / (1 / r2 * temp + 1 / r1 * (1 - temp));
        }
        else if(ff_vf[cf]<0)
        {
                Scal temp=1e-6;
                t.ff_resist[cf] = 1 / (1 / r2 * temp + 1 / r1 * (1 - temp));
        }
        else
        {
                t.ff_resist[cf] = 1 / (1 / r2 * ff_vf[cf] + 1 / r1 * (1 - ff_vf[cf]));
        }
```

2. Compile the **Aphros** code by following the instructions described [here](https://github.com/cselab/aphros).

## Run

Navigate to the `setup` folder and run with:

```
make cleanrun
```

The simulation parameters including *bubble nucleation radius*, *flow rate*, and *current* can be modified in `setup` file.

In `std.conf`, the following parameter can be modified to simulate with or without bubble coalescence.

With bubble coalescence:

```
set string advection_solver vof
```

Without bubble coalescence:

```
set string advection_solver vofm
```

