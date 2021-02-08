Buoyancy
========

In aquatic oil spills, we can treat the oil slick as being composed of small droplets of oil suspended
in the surrounding water. These droplets are subjected to buoyancy, the upward force
imposed on the droplet by the water. Buoyancy is greater deeper in the water column.

Terminology
-----------



Terminal Velocity Equations
---------------------------

As discussed in :doc:`droplet_shapes`, oil particles are spherical at small sizes, ellipsoid at
intermediate sizes, and form a spherical cap at larger sizes. Each size range requires
a variation on the terminal velocity equation. The following equations come from Clift (1978).

Spherical Small Size (diameter less than 1 mm)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. math::

   U_T = \frac{R \mu}{\rho d}

Where the Reynold's number :math:`R` is calculated from the Best number :math:`N_D` as follows:

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



References
----------

Clift, R., J. R. Grace, and M. E. Weber, 2005, `Bubbles, Drops, and Particles <https://www.google.com/books/edition/_/UUrOmD8niUQC?hl=en>`_.

Zheng Li, and Yapa Poojitha D., 2000, `Buoyant Velocity of Spherical and Nonspherical
Bubbles/Droplets <https://ascelibrary.org/doi/abs/10.1061/(ASCE)0733-9429(2000)126:11(852)>`_:
Journal of Hydraulic Engineering, 126, 852–854.
