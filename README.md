# TB Transmission in Kenya: What Will It Take to Meet the 2030 Target?

Kenya has one of the highest TB burdens in Africa. Despite decades of investment in its national TB programme, incidence remains far above the WHO End TB 2030 target. This project asks a simple but consequential question: under realistic intervention scenarios, how close can Kenya actually get?

To answer it, I built and calibrated a compartmental transmission model using WHO Kenya incidence data and Kenya-specific parameters from the published literature. The model is deliberately simple with five compartments, fixed parameters, single-point calibration. I think that honesty about what a model can and cannot do is more useful than complexity that obscures 
its limitations.

## The Questions This Project Addresses

1. Is Kenya's TB epidemic self-sustaining without intervention, and what does the basic reproduction number tell us about that?

2. What does Kenya's observed post-2000 incidence decline tell us about what actually drove it, biology or programme improvement?

3. Which programmatic parameter: treatment initiation, transmission reduction, or both combined, offers the greatest leverage for reducing TB burden by 2030?

4. Which biological and programmatic factors most strongly determine TB incidence in Kenya, and what does that mean for where the health system should focus?

## What I Found

Kenya's calibrated R0 is 1.10. TB persists without intervention. The model cannot reproduce the observed post-2000 decline in incidence, and I think that is informative rather than a failure. It suggests the decline was driven by improving treatment coverage over time, not by changes in transmission biology. A fixed-parameter model is not designed to capture that.

None of the three intervention scenarios reaches the End TB 2030 target of 76 per 100,000. Even combined treatment scale-up and transmission reduction projects 216 per 100,000 by 2030, nearly three times the target. The sensitivity analysis points to why: treatment initiation (PRCC = -0.61) is the strongest lever available, but pulling it harder 
is not enough on its own.

| Scenario | 2030 Projection |
|---|---|
| Status quo | 412 per 100,000 |
| Treatment scale-up (tau +30%) | 299 per 100,000 |
| Combined intervention | 216 per 100,000 |
| End TB target | 76 per 100,000 |

## How the Model Works

Five compartments: Susceptible, Latent, Infectious, On Treatment, Recovered. I modelled treatment explicitly rather than folding it into a generic recovery rate. This lets me calibrate directly to Kenya's NTLD programme indicators and run treatment-specific scenarios.

Calibration uses a 50-year burn-in to reach endemic equilibrium before fitting the transmission rate to the 2000 baseline incidence. Scenarios are grounded in Kenya's NTLD Strategic Plan targets, not arbitrary parameter changes. Uncertainty bounds on projections reflect parameter uncertainty in the transmission rate.

## Repository Structure
```
tb-kenya-transmission-model/
├── R/
│   ├── model.R          # ODE system and R0 function
│   ├── calibration.R    # Beta calibration to WHO Kenya data
│   └── sensitivity.R    # PRCC sensitivity analysis
├── data/
│   ├── who_kenya_tb.csv       # WHO Global TB Report
│   └── kenya_tb_who.csv       # Kenya incidence 2000-2022
├── outputs/
│   └── tb_kenya_model.html    # Rendered Quarto report
└── tb_kenya_model.qmd         # Main analysis document
```
## Data Sources

- WHO Global TB Report (who.int/tb/data)
- Kenya NTLD Programme Annual Reports 2021 and 2022
- KNBS Statistical Abstract 2023
- Kirimi et al. (2024)- model parameters

## Limitations Worth Knowing

Fixed parameters across the full study period. Single latent compartment. 
No HIV co-infection. Single-point calibration. All acknowledged in the 
report with specific suggestions for future work.

## Tools

R, deSolve, tidyverse, Quarto, Git

