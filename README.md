# GEEmix
  Public Available Code and Data to Reproduce Analyses in "Marginal Analysis of Bivariate Mixed Responses with Measurement Error and Misclassification."

## File Structure

### Data Analysis

* Setting 1: Complicated measurement error process
  * [GEEmix](https://github.com/QihuangZhang/GEEmix/blob/main/code/DataAnalysis/DataAnalysis_905_IVoptimcompcomp.R)
  * [GMM](https://github.com/QihuangZhang/GEEmix/blob/main/code/DataAnalysis/DataAnalysis_905_IVcompGMM.R)
* Setting 2: Simplified measurement error process
  * [GEEmix](https://github.com/QihuangZhang/GEEmix/blob/main/code/DataAnalysis/DataAnalysis_905_IVoptimcompsimp.R)
  * [GMM](https://github.com/QihuangZhang/GEEmix/blob/main/code/DataAnalysis/DataAnalysis_905_IVsimpGMM.R)

### Simulation

* Simulation 1: mismeasurement parameters are known
  * Case 1 - linear measurement error process: [Simulation_INS2_2](https://github.com/QihuangZhang/GEEmix/blob/main/code/Simulation/Simulation_INS2_2.R)
  * Case 2 - non-linear measurement error process: [Simulation_INS12_2](https://github.com/QihuangZhang/GEEmix/blob/main/code/Simulation/Simulation_INS12_2.R)
* Simulation 2: validation data are available
  * Internal Validation: [Simulation_INS3](https://github.com/QihuangZhang/GEEmix/blob/main/code/Simulation/Simulation_INS3.R)
  * External Validation: [Simulation_INS4](https://github.com/QihuangZhang/GEEmix/blob/main/code/Simulation/Simulation_INS4.R)
  * No measurement error in response: [Simulation_INS8](https://github.com/QihuangZhang/GEEmix/blob/main/code/Simulation/Simulation_INS8.R)
* Simulation 3: relative efficiency of the proposed estimators with internal validation data
  * Weighted (Internal) Estimator: [Simulation_INS5](https://github.com/QihuangZhang/GEEmix/blob/main/code/Simulation/Simulation_INS5.R)
  * Generalized Method of Moment (GMM): [Simulation_INS10_2](https://github.com/QihuangZhang/GEEmix/blob/main/code/Simulation/Simulation_INS10_2.R)
  * External validation (in comparison): [Simulation_INS6](https://github.com/QihuangZhang/GEEmix/blob/main/code/Simulation/Simulation_INS6.R)
  * Figure 2 and 3: [Simulation_INS7](https://github.com/QihuangZhang/GEEmix/blob/main/code/Simulation/Simulation_INS7.R)