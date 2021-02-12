Lagrangian Particle Tracking
============================

Lagrangian particle tracking (LPT) simulates the flow of an oil slick in water by modeling the slick as a collection of discrete droplets. 

Velocity
--------

2D Surface Velocity
+++++++++++++++++++

The velocity of a surface droplet is influenced both by current velocity in the water and wind velocity. We can write
the droplet velocity, :math:`\boldsymbol{u_d}`, at a position :math:`\boldsymbol{x_d} = (x_d,\,y_d)` as a combination of 
current velocity, :math:`\boldsymbol{u_c}`, and wind velocity, :math:`\boldsymbol{u_w}`:

.. math::
   \boldsymbol{u_d}(\boldsymbol{x_d}) = F_c\boldsymbol{u_c}(\boldsymbol{x_d}) + F_w(z_d) \boldsymbol{u_w}(\boldsymbol{x_d})

where :math:`F_c` and :math:`F_w` are the current and wind modifiers, respectively. 

3D Velocity
+++++++++++

When extending to 3D, we incorporate :ref:`buoyancy <buoyancy>` into the :math:`z_d` component of :math:`\boldsymbol{x_d} = (x_d,\,y_d,\,z_d)`:

.. math::
   u_{d_z}(\boldsymbol{x_d}) = F_c u_{c_z}(\boldsymbol{x_d}) + F_w(z_d) u_{w_z}(\boldsymbol{x_d}) + u_T(\boldsymbol{x_d})

where :math:`u_{c_z}` and :math:`u_{w_z}` are the :math:`z` components of :math:`\boldsymbol{u_d}` and :math:`\boldsymbol{u_w}`, respectively, and 
:math:`u_T` is the :ref:`bouyant velocity <buoyancy>`.


Diffusion
---------

References
----------

Dietrich, J. C., C. J. Trahan, M. T. Howard, J. G. Fleming, R. J. Weaver, S. Tanaka, L. Yu, R. A. Luettich, C. N. Dawson, 
J. J. Westerink, G. Wells, A. Lu, K. Vega, A. Kubach, K. M. Dresback, R. L. Kolar, C. Kaiser, and R. R. Twilley, 2012, 
`Surface trajectories of oil transport along the Northern Coastline of the Gulf of Mexico <https://www.sciencedirect.com/science/article/abs/pii/S0278434312000799>`_: 
Continental Shelf Research, 41, 17–47.

Wang, S.-D., Y.-M. Shen, Y.-K. Guo, and J. Tang, 2008, 
`Three-dimensional numerical simulation for transport of oil spills in seas <https://www.sciencedirect.com/science/article/pii/S0029801807002582>`_: 
Ocean Engineering, 35, 503–510.

