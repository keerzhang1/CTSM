# Software Design Documentation

<span style="font-size:larger;"><b>Project Name</b></span>
<!-- Note you can add this to an issue or pull request in github. That is the normal usage of this template. -->

**Date:**: Updated on June 28, 2025

**Written By**: Keer Zhang

## Introduction
---------------------------------------
*What is the goal of the software? What is the problem statement?*

If there are any specification documents, link them in Appendix.

We added urban tree in the CLM urban model. 

## Solutions
---------------------------------------
*Section should include alternative implementations/solutions*

Is it feasible? How much effort does it need for each approach? Pros/cons of each approach.

Document alternatives, why you made the decision and how it will affect the team and project.

Tree layers: consider two (above-roof and below-roof) tree layers 
It's easier to only implement one tree layer, but we'd like to consider more complex energy 
distribution at the lower and upper part of tree canopy. It's important to consider tree shading 
above roofs.

Radiative process:
1. We use Monte Carlo Ray Tracing to calculate view factors from one surface to another

## Design Considerations
---------------------------------------
Describe the issues that need to be addressed before creating a design solution.



### Assumptions and Dependencies:
Describe any assumptions that may be wrong or any dependencies on other things



### General Contraints:
Describe any constraints that could have an impact on the design of the software.


## Design and Architecture
---------------------------------------

### System diagram or flowchart
Interaction diagram of various inputs, outputs, sub systems and dependencies.


### Algorithm or Pseudo code for main components
Describe your logic in this section

## Rollout Plan
---------------------------------------
Define the roll-out phases and tests you plan to do

## Appendix
---------------------------------------
References, links to additional documentation

## Review Sign-off
---------------------------------------
* Reviewer(s):

*Sign-off Completed on YYYY-MM-DD*
