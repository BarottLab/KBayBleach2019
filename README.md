# KBayBleach2019
This repository contains the data and code to accompany the manuscript: 

Innis, T., Allen-Waller, L., Brown, K.T., Sparagon, W., Carlson, C., Kruse, E., Huffmyer, A., Nelson, C., Putnam, H., and Barott, K. (2021). Marine heatwaves depress metabolic activity and impair cellular acid-base homeostasis in reef-building corals regardless of bleaching susceptibility

Abstract:
Ocean warming is causing global coral bleaching events to increase in frequency, resulting in widespread coral mortality and disrupting the function of coral reef ecosystems. However, even during mass bleaching events, many corals resist bleaching despite exposure to abnormally high temperatures. While the physiological harm of bleaching has been well documented, the consequences of heat stress for bleaching resistant individuals are not well understood. In addition, much remains to be learned about how heat stress affects cellular level processes that may be overlooked at the organismal level, yet are crucial for coral performance in the short term and ecological success over the long term. Here we compared the physiological and cellular responses of bleaching resistant and bleaching susceptible corals throughout the 2019 marine heatwave in Hawai‘i, a repeat bleaching event that occurred four years after the previous regional event. Relative bleaching susceptibility within species was consistent between the two bleaching events, yet corals of both resistant and susceptible phenotypes exhibited pronounced metabolic depression during the heatwave. At the cellular level, bleaching susceptible corals had lower intracellular pH than bleaching resistant corals at the peak of bleaching for both symbiont-hosting and symbiont-free cells, indicating that disruption of acid-base homeostasis was worse in bleaching susceptible individuals. Notably, cells from both phenotypes were unable to compensate for experimentally induced cellular acidosis, indicating that acid-base regulation was significantly impaired at the cellular level even in bleaching resistant corals and even in cells containing symbionts. Thermal disturbances may thus have substantial ecological consequences, as even small reallocations in energy budgets to maintain homeostasis during stress can negatively affect fitness. These results suggest concern is warranted for corals coping with ocean acidification alongside ocean warming, as the feedback between temperature stress and acid-base regulation may further exacerbate the physiological effects of climate change.

**Repository contents:**

**data/:** Contains environmental and field collection data, physiology measurements and acidosis experiment measurements. 
* raw phys data/: raw measurements for all physiology variables
* pHi/: acidosis experiment measurements
* benthos/: visual bleaching response
* temp/: temperature data obtained from NOAA sensor on Moku o loe and DHW/BI summaries
* PhysData.csv: physiology measurement summary file

**analysis/:**
* pHi_analysis_2019KBay.Rmd: acidosis experiment analysis
* pi_curves.Rmd: generating PI curves from PreSens output
* DHW_BI.Rmd: calculations for reef-wide bleaching index and DHW
* mixedmodels.Rmd: statistical analysis for physiology measurements
* pi_models.Rmd: PI curve model analysis
