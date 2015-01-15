A Brief History of Geoscience Australia's Image Processor
=========================================================

The gaip code base has evolved through 4 major iterations evolving from a system starting as a system that was essentially a set of shell scripts that ran Fortran 77 executables to a job runner framework that utilises Luigi.  As part of the 4th iteration, a large portion of the Fortran 77 and C code has been re-written in Python, to take advantage of the larger codebase available through NumPy, SciPy, and GDAL.
Original implementations came about as part of the Unlocking the Landsat Archive (ULA) project.

ULA1
----

A description of the original system and some notes on its implementation can be found in TRIM 2009/4427. Of particular interest, D2010-36176 :download:`Landsat NBAR Product Generation System <auxiliary/Landsat NBAR Product Generation System.docx>` which gives a broad overview of the system and references D2010-16918 :download:`NBA Project Approach document <auxiliary/NBA Project Approach document.docx>`  which describes some of the business needs and the structure of the project team (members of the current ULA team note further people whom were involved in the project and are noted in the list below).

From :download:`auxiliary/Landsat NBAR Product Generation System <Landsat NBAR Product Generation System.docx>`:

*The key scientific algorithm of the system was implemented in FORTRAN by Fuqin Li while the main driver and controlling modules were in Python. Nine Shell scripts, developed by Fuqin Li and Lan-Wei Wang, were used to communicate between the Python and FORTRAN modules. A couple of C programs were developed for the extraction of ancillary parameters and the generation of quick-look JPG images. The system can be run on a stand-alone Linux environment and has been successfully tested on acr111, ausimages4 and ANU NCI Linux machines.*

**The people involved in ULA1 were:**

* `Fei Zhang <mailto:fei.zhang@ga.gov.au>`_ (ULA team leader),
* `Frank Fu <mailto:frank.fu@ga.gov.au>`_ (developer),
* `Paul Gardner <mailto:paul.gardner@ga.gov.au>`_ (developer),
* `Lan-Wei Wang <mailto:lan-wei.wang@ga.gov.au>`_ (developer and business development),
* `Fuqin Li <mailto:fuqin.li@ga.gov.au>`_ (scientist),
* `Leo Lymburner <mailto:leo.lymburner@ga.gov.au>`_ (scientist),
* `Medhavy Thankappan <medhavy.thankappan@ga.gov.au>`_ (project manager), and
* `Wnjun Wu <wenjun.wu@ga.gov.au>`_ (senior supplier).

Note that not all these people worth within the ULA team or NEO any more.

ULA2
----

The second iteration of ULA development largely comprised of numerical modifications that generalised the calculations to non rectangular regions.

The main person responsible involved in this iteration was `Roger Edburg <mailto:roger.edburg@ga.gov.au>`_.

**The people involved in ULA2 were:**

* `Simon Oliver <mailto:simon.oliver@ga.gov.au>`_ (ULA team leader),
* `Roger Edburg <mailto:roger.edburg@ga.gov.au>`_, and
* `Paul Gardner <mailto:paul.gardner@ga.gov.au>`_.

ULA3
----

The ULA image processor implements a re-architected version of the original ULA NBAR processor within an internal workflow management framework which manages the tasks and sub-tasks in a structured manner. The framework was designed to be able to incorporate other operations (such as PQA and FC) at a later date.

The new system uses exactly the same algorithms as the previous ULA NBAR processor, but has been refactored to provide a more logical and efficient structure. The new system also implements multi-threading which can significantly reduce scene turnaround time on machines with sufficient resources.The NBAR processing algorithm is constantly evolving, so it was necessary to develop a system which would permit new tasks to be easily integrated while minimising the programming overhead. Re-arrangement of the workflow can be performed by editing configuration files rather than Python code.

**The people involved with ULA3 were:**

* `Matthew Purss <mailto:matthew.purss@ga.gov.au>`_ (ULA team leader),
* `Alex Ip <mailto:alex.ip@ga.gov.au>`_ (developer),
* `Roger Edburg <mailto:roger.edburg@ga.gov.au>`_ (developer),
* `Paul Gardner <mailto:paul.gardner@ga.gov.au>`_ (developer),
* `Joshua Sixsmith <mailto:joshua.sixsmith@ga.gov.au>`_ (scientist/developer),
* `Fuqin Li <mailto:fuqin.li@ga.gov.au>`_ (scientist), and
* `Simon Knapp <mailto:simon.knapp@ga.gov.au>`_.

gaip
----

The gaip module was re--architected to make use of F2Py for some of the Fortran 77 routines, make use of existing 3rd party libraries, reduce the overall memory footprint, modularise the code to make a simpler interface and encorages portability, as well as simplify the workflow using Luigi as the framework. Positive outcomes was a more simplified and maintainable codebase.  The Pixel Quality (PQ) and Fractional Cover (FC) modules were refactored into the Luigi framework as well as reductions in overall memory use.
The NBAR framework now outputs 3 different reflectance products:
1. Lambertian
2. BRDF corrected
3. BRDF and Terrain corrected
The BRDF correction also differs from the ULA3 BRDF correction, in that a more accurate model is used.

**The people involved with gaip were:**

* `Lan-Wei Wang <mailto:lan-wei.wang@ga.gov.au>`_ (Stream Lead),
* `Dale Roberts <mailto:dale.roberts@ga.gov.au>`_ (Developer),
* `Steven Ring <mailto:steven.ring@ga.gov.au>`_ (Developer),
* `Josh Sixsmith <mailto:joshua.sixsmith@ga.gov.au>`_ (Developer), and
* `Fuqin Li <mailto:fuqin.li@ga.gov.au>`_ (Scientist).
