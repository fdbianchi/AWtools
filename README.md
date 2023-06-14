# AWtools

A toolbox aimed to design anti-windup (AW) strategies for MIMO plants. 

The design tools are based on the traditional two-step approaches, in which the AW compensator is designed after a controller is tuned without considering saturation. 

AWtools covers LTI and LPV plants and is based on the comprime factor formulation introduced in:
- M. C. Turner and I. Postlethwaite, ‘A new perspective on static and low order anti-windup synthesis’, International Journal of Control, vol. 77, no. 1, pp. 27–44, 2004, doi: 10.1080/00207170310001640116.
- S. Skogestad and I. Postlethwaite, Multivariable feedback control: analysis and design, 2. ed. Chichester: Wiley, 2005.

Examples of how to use it can be found in the folder **demo**.

The **awtools.pdf** provides a detail description of the main function **awsyn**

