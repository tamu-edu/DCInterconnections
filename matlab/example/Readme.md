This repository provides a **ready-to-run demonstration** for simulating the grid impact of a **1 GW data center load**.  

To keep the example lightweight and easy to reproduce, the simulation scope has been intentionally simplified compared to the full-scale study.

The code is intended to showcase an **end-to-end workflow** (reliability screening → price/congestion metrics), rather than to represent the full scale simulation.

---
## Simulation Configuration
- **Target buses:** Bus 30–50  
- **Data center size:** 1 GW total load  
- The 1 GW data center demand is added to each candidate target bus for evaluation.

---

## How to Run

### Step 1 — run reliability screening to obtain the qualified nodes
Run:
```matlab
Reliability_screen.m
```



### Step 2 — run price/congestion analysis to obtain the price-congestion metrics
Run:
```matlab
Flex_DC_publish.m
```
