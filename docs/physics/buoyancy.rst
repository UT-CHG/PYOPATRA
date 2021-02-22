.. _buoyancy:

Buoyancy
========

In aquatic oil spills, we can treat the oil slick as being composed of small droplets of oil suspended
in the surrounding water. These droplets are subjected to buoyancy, the upward force
imposed on the droplet by the water. Buoyancy is greater deeper in the water column.

Terminology
-----------

We say that a *droplet* (e.g., of oil) is contained within a *continuous fluid* (e.g., water).

Buoyancy for particles involves several physical quantities:

- :math:`U`, the velocity of a droplet, and :math:`U_T`, the terminal velocity of a droplet.
- :math:`\nu`, the kinematic viscosity of the continuous fluid within which the 
  droplets are contained.
- :math:`d`, the diameter of the droplet, for spheres, or :math:`d_e`, the effective 
  diameter for ellipsoids or spherical caps (i.e., the diameter of an equivalent volume sphere).
- :math:`R`, Reynolds Number: the ratio of inertial forces to viscous forces. 
  typically :math:`Ud/\nu` for particles.
- :math:`\rho`, the density of the continuous fluid, and :math:`\Delta \rho`,
  the difference between the density of the continuous fluid and the 
  density of the droplet.
- :math:`g`, the acceleration due to gravity. :math:`9.8 \, \text{m/s}^2` in this package.
- :math:`\mu`, the viscosity of the continuous fluid, and :math:`\mu_w` the viscosity of water
  from Braida (0.0009 kg/ms).
- :math:`M`, Morton number, and :math:`E_O`, Eötvös number, used to characterize the shape of ellipsoid 
  spherical cap particles.
- :math:`\sigma`, the interfacial tension.


Water Viscosity
---------------

Water viscosity plays an important role in buoyancy. However, water viscosity is
dependent on density, :math:`\rho`, and temperature, :math:`T`, both of which depend
on depth in bodies of water. So, to calculate water viscosity, we must take into account
the depth of the particle in question.

Water Viscosity With Respect to Density
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Ocean water contains salt, which can impact density as well as temperature.
We follow Huber et al (2009) to calculate water viscosity. In essence, the dimensionless
water velocity is calculated as

.. math::
   \bar{\mu} = \bar{\mu}_0(\bar{T}) \times \bar{\mu}_1(\bar{T},\,\bar{\rho}))

where :math:`\bar{\mu}_0` is viscosity when density is 0 (in the limit) and :math:`\bar{\mu}_1` is viscosity due to increasing
density. Furthermore, :math:`\bar{T}` and :math:`\bar{\rho}`
represent dimensionless quantities:

.. math::
   \begin{align}
      \bar{T} &= \frac{T}{T^\star} = \frac{T}{647.096 \text{ K}} \\
      \bar{\rho} &= \frac{\rho}{\rho^\star} = \frac{\rho}{322.0 \text{ kg m}^{-3}}
   \end{align}

Viscosity at the Zero Density Limit
+++++++++++++++++++++++++++++++++++

We use the following equation from Huber et al (2009):

.. math::
   \bar{\mu}_0 = \frac{100 \, \bar{T}}{1.67752 + \frac{2.20462}{\bar{T}} + \frac{0.6366564}{\bar{T}^2} - \frac{0.241605}{\bar{T}^3}}.

The coefficients were obtained by orthogonal distance regression on observational data.

Viscosity due to Increasing Density
+++++++++++++++++++++++++++++++++++

This term depends on both dimensionless temperature :math:`\bar{T}`, and on dimensionless density,
:math:`\bar{\rho}`:

.. math::
   \bar{\mu}_1(\bar{T},\,\bar{\rho}) = \exp\left[ \sum^5_{i=0}\left( \frac{1}{\bar{T}} - 1 \right)^i \sum^6_{j=0} H_{ij} \left( \bar{\rho} - 1 \right)^j \right]

where :math:`H_{ij}` are determined by the following table:

.. list-table:: Viscosity Coefficients
   :header-rows: 1

   * - :math:`i`
     - :math:`j`
     - :math:`H_{ij}`
   * - 0
     - 0
     - 0.520094
   * - 1
     - 0
     - 0.0850895
   * - 2
     - 0
     - -1.08374
   * - 3
     - 0
     - -0.289555
   * - 0
     - 1
     - 0.222531
   * - 1
     - 1
     - 0.999115
   * - 2
     - 1
     - 1.88797
   * - 3
     - 1
     - 1.26613
   * - 5
     - 1
     - 0.120573
   * - 0
     - 2
     - -0.281378
   * - 1
     - 2
     - -0.906851
   * - 2
     - 2
     - -0.772479
   * - 3
     - 2
     - -0.489837
   * - 4
     - 2
     - -0.257040
   * - 0
     - 3
     - 0.161913
   * - 1
     - 3
     - 0.257399
   * - 0
     - 4
     - -0.0325372
   * - 3
     - 4
     - 0.0698452
   * - 4
     - 5
     - 0.00872102
   * - 3
     - 6
     - -0.00435673
   * - 5
     - 6
     - -0.000593264

Values not in the table are 0.0.

Pure Water Viscosity with Respect to Temperature
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For pure water we make a simplifying assumption and assume we are at atmospheric pressure. Though pressure does increase
with depth in the ocean, the contribution of the increased pressure to density is small. We use the following equation
from Patek et al (2009):

.. math::
   \bar{\mu} = 280.68 \bar{T}^{-1.9} + 511.45 \bar{T}^-{7.7} + 61.131 \bar{T}^{-19.6} + 0.45903 \bar{T}^{-40.0}

where :math:`\bar{T} = T / 300` K. This equation is only recommended in the temperature range of 253.15 K to 383.15 K
(-20 C to 110 C).


Terminal Velocity Equations
---------------------------

As discussed in :doc:`droplet_shapes`, oil particles are spherical at small sizes, ellipsoid at
intermediate sizes, and form a spherical cap at larger sizes. Each size range requires
a variation on the terminal velocity equation. The following equations come from Clift (1978).

Spherical Small Size (diameter less than 1 mm)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. math::

   U_T = \frac{R \mu}{\rho d}

Where the Reynold's number, :math:`R`, is calculated from the Best number, :math:`N_D`, as follows:

.. math::

   N_D = \frac{4\rho \Delta \rho g d^3}{3 \mu^2}

.. list-table:: Reynolds Number Calculation
   :header-rows: 1

   * - Range
     - Correlation
   * - :math:`N_D \leq 73`
     - :math:`N_D / 24 - 1.7569 \times 10^{-4} N_D^2 + 6.9252 \times 10^{-7} N_D^3 - 2.3027 \times 10^{-10} N_D^4`
   * - :math:`73 < N_D \leq 580`
     - :math:`\log R = -1.7095 + 1.33438 \log N_D - 0.11591 (\log N_D)^2`
   * - :math:`580 < N_D \leq 1.55 \times 10^7`
     - :math:`\log R = -1.81391 + 1.34671 \log N_D - 0.12427 (\log N_D)^2 + 0.006344 (\log N_D)^3`


Ellipsoid (Intermediate Size) and Spherical Caps (Large Size)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We require the Morton number and Eötvös numbers:

.. math::
   
   M = \frac{g \mu^4 \Delta \rho}{\rho^2 \sigma^3}

   E_O = \frac{g \Delta \rho d^2_e}{\sigma}.

Ellipsoid (:math:`M < 10^{-3}` and :math:`E_O \leq 40`)
+++++++++++++++++++++++++++++++++++++++++++++++++++++++

We have:

.. math:: 
   U_T = \frac{\mu}{\rho d_e} M^{-0.149} (J - 0.857).
   :label: ellipsoid
   
   

We calculate an intermediate value, :math:`H`:

.. math::

   H = \frac{4}{3} E_O M^{-0.149}(\mu / \mu_w )^{-0.14}.

For :math:`2 < H \leq 59.3`,

.. math::

   J = 0.94 H^{0.757}.

For :math:`H > 59.3`,

.. math::

   J = 3.42 H^{0.441}.


Spherical Cap (:math:`E_O > 40`)
++++++++++++++++++++++++++++++++

In this range, we simply have

.. math:: 
   U_T = 0.711 \sqrt{g d_e \Delta \rho / \rho}
   :label: spherical-cap

   
   


Critical Diameter
+++++++++++++++++

Zheng and Yapa (2000) introduce an approximate method to find the critical diameter which separates the ellipsoid droplets from spherical caps.
They assume that Equations :eq:`ellipsoid` and :eq:`spherical-cap` represent straight lines in logarithmic coordinates, and then locate the intersection
of the two lines to find the critical diameter:

.. math::
   \log d_c = \frac{b_2 - b_1}{a_1 - a_2}

where, for points :math:`(x_1,\,y_1)` and :math:`(x_2,\,y_2)` of :math:`(\log d_e,\,\log U_T)` when :math:`H=59.3` and :math:`d_e = 15` mm, respectively,

.. math::
   \begin{align}
      a_1 &= 0.5 \\
      b_1 &= \log \left[0.711 (g \Delta \rho / \rho)^{1/2}\right] \\
      a_2 &= \frac{y_2 - y_1}{x_2 - x_1} \\
      b_2 &= y_1 - a_2 x_1.
   \end{align}



References
----------

Clift, R., J. R. Grace, and M. E. Weber, 2005, `Bubbles, Drops, and Particles <https://www.google.com/books/edition/_/UUrOmD8niUQC?hl=en>`_.

Huber, M. L., R. A. Perkins, A. Laesecke, D. G. Friend, J. V. Sengers, M. J. Assael, I. N. Metaxa, E. Vogel, R. Mareš, and K. Miyagawa, 2009, `New international formulation for the viscosity of H2O <https://aip.scitation.org/doi/full/10.1063/1.3088050>`_: Journal of Physical and Chemical Reference Data, 38, 101–125.

Zheng Li, and Yapa Poojitha D., 2000, `Buoyant Velocity of Spherical and Nonspherical
Bubbles/Droplets <https://ascelibrary.org/doi/abs/10.1061/(ASCE)0733-9429(2000)126:11(852)>`_:
Journal of Hydraulic Engineering, 126, 852–854.
