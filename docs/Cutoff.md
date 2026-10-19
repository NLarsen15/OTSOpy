# Cutoff

![Cutoff Diagram](assets/cutoffplot.png "Computation of the Oulu neutron monitor effective cut-off rigidity using the IGRF 2000 epoch and TSY89 model with kp index = 0. Penumbra is shown by the forbidden and allowed trajectories being black and white, respectively. The upper and lower cut-off values (Ru and Rl) are denoted in the legend, from which the effective cut-off (Rc) is computed.")

<figure markdown>
<video autoplay loop muted playsinline controls poster="../assets/OTSO_oulu_penumbra_poster.jpg" style="width: 100%;" title="Oulu neutron monitor penumbra traced with OTSO">
  <source src="../assets/OTSO_oulu_penumbra.webm" type="video/webm">
  <source src="../assets/OTSO_oulu_penumbra.mp4" type="video/mp4">
</video>
<figcaption style="max-width: 100%;">Vertical cosmic rays backtraced from the Oulu neutron monitor at decreasing rigidity using the IGRF and TSY89 models with kp = 0 for 1 January 2024 12:00 UT (the OTSO defaults). These conditions differ from the figure above (IGRF 2000 epoch), so the penumbra structure and cut-off values are not the same, mainly because the geomagnetic field has changed between the two epochs. Above the penumbra all trajectories are allowed (blue); within it allowed and forbidden (orange) bands alternate and the trajectories become increasingly complex; below it all are forbidden. The bar fills in the allowed and forbidden pattern at 0.005 GV resolution, and the upper, effective and lower cut-off rigidities (Ru, Rc and Rl) are marked at the end.</figcaption>
</figure>

::: OTSO.cutoff
    options:
      show_source: false
