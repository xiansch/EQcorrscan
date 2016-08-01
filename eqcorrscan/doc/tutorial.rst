EQcorrscan tutorial
===================
Welcome to EQcorrscan - this package is designed to compute earthquake detections
<<<<<<< HEAD
using a paralleled matched-filter network cross-correlation routine.  The inner
loop of this package is the cross-correlation of templates of seismic data
with day-long seismic data.  This inner function is the openCV.match_template
function - this appears to be a well optimized cross-correlation function, and
is written in c++.  Cross-correlations are computed in the frequency domain
for large datasets, for which a day of seismic data usually qualifies.

Before continuing with this tutorial please check that you have installed all
the pre-requisite modules, as not all will be installed by the setup.py file.
The list of these is in the Introduction section of this documentation.

As you will see, this package is divided into two main sub-modules, the
Core and Utils sub-modules.  The Core sub-module contains the main, high-level
functions:
=======
using a paralleled matched-filter network cross-correlation routine, and analyse the
results.

Before continuing with this tutorial please check that you have installed all
the pre-requisite modules, as not all will be installed by the setup.py file.
The list of these is in the :doc:`Introduction <intro>` section of this documentation.

As you will see, this package is divided into two main sub-modules, the
:doc:`core </core>` and :doc:`utils </utils>` sub-modules.
The core sub-module contains the main, high-level functions:
>>>>>>> upstream/master

:bright_lights:
        A brightness based template detection routine;
:template_gen:
        A series of routines to generate templates for match-filter detection
<<<<<<< HEAD
        from continuous or cut data, with pick-times defined either manually, or from a
        *Seisan* s-file;
:match_filter:
        The main matched-filter routines, this is split into several
        smaller functions to allow python based parallelisation;
:lag_calc:
        Routines for calculating optimal lag-times for events detected
        by the match-filter routine, these lags can then be used to define new picks
        for high accuracy relocations. *Under-development*

The Utils sub-module contains useful, but small functions.  These functions are
rarely cpu intensive, but perform vital operations, such as reading *Seisan* s-files,
finding peaks in noisy data, converting a seisan database to hypoDD formatted
files and computing cross-correlations between detections for hypoDD (a double
difference relocation software), calculating magnitudes, clustering detections,
stacking detections, making pretty plots, and processing seismic data in the
same way repeatedly using *Obspy*'s functionality.

What follows is a (hopefully) expanding set of tutorials that should take you
=======
        from continuous or cut data, with pick-times either defined manually,
        or defined in event files;
:match_filter:
        The main matched-filter routines, this is split into several
        smaller functions to allow python-based parallel-processing;
:lag_calc:
        Routines for calculating optimal lag-times for events detected
        by the match-filter routine, these lags can then be used to define new picks
        for high accuracy re-locations. *Under-development*

The :doc:`utils </utils>` sub-module contains useful, but small functions.
These functions are rarely cpu intensive, but perform vital operations, such
as reading *Seisan* s-files (:doc:`sfile_util </submodules/utils.sfile_util>`),
finding peaks in noisy data (:doc:`findpeaks </submodules/utils.findpeaks>`),
converting a seisan database to hypoDD formatted files and computing cross-correlations between
detections for hypoDD (a double difference relocation software)
(:doc:`catalog_to_dd </submodules/utils.catalog_to_dd>`), calculating
magnitudes (:doc:`mag_calc </submodules/utils.mag_calc>`),
clustering detections (:doc:`clustering </submodules/utils.clustering>`),
stacking detections (:doc:`stacking </submodules/utils.stacking>`),
making pretty plots (:doc:`plotting </submodules/utils.plotting>`),
and processing seismic data in the same way repeatedly using *Obspy*'s
functionality (:doc:`pre_processing </submodules/utils.pre_processing>`).

What follows is an expanding set of tutorials that should take you
>>>>>>> upstream/master
through some of the key functionality of the EQcorrscan package.

.. toctree::
  :numbered:
  :titlesonly:

  tutorials/template-creation.rst
  tutorials/matched-filter.rst
<<<<<<< HEAD
  tutorials/lag-calc.rst
=======
>>>>>>> upstream/master
  tutorials/mag-calc.rst
  tutorials/clustering.rst
