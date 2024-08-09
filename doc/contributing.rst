Contributing
==============

In this section, we lay out rules and recommendations for contributing to this project. A set of rules (a coding style guide) is essential for any development team to ensure code readability and integrity. As a developer, you should pause to go through either all of the rules or ones that are relevant to you depending on what you are currently working on, so that when you commit new code it will pass all linter, regression, and other tests.

It is worth mentioning that all guidelines listed in section are quite standard and used by the majority of the Python and C++ communities. Adopting and getting used to these will hopefully also benefit you with any future projects.


Git
---

This section is a very minimal tutorial on git. It is strictly speaking not part of a coding style but aims to minimize the risk of making bad commits or losing data. Do pay attention to the paragraph on commit messages though -- following the correct style will make the commit history much more readable!


Basic usage
^^^^^^^^^^^^

Navigate into any source directory of the project and run

.. code-block:: bash

   git status

to see which files have been modified. Run

.. code-block:: bash

   git add

on each file you want to commit. Similarly, run

.. code-block:: bash

   git rm

on each file you want to remove from version control (don’t remove them with ``rm``). In order to commit, issue

.. code-block:: bash

   git commit

which prompts you with the commit message before the actual commit is performed. Attention: never run ``git commit -a`` unless you are an experienced Git user!

The basics of how to write commit messages are well explained in this blog post: https://chris.beams.io/posts/git-commit. In short, start with a summary line consisting of no more than 50 characters, not followed by a period. Leave a blank line followed by further description if necessary. For small commits, just the summary line may be sufficient. Write the whole commit message in the imperative tense (i.e. "Fix typo" not "Fixed typo"). These might seem like arbitrary rules but getting in the habit of creating quality commit messages makes using and collaborating with Git a lot easier! For instance, go ahead and issue ``git log --oneline`` in the source directory. You will notice that only the header line is taken from the commit messages. If these are well written it makes searching for a particular commit a lot easier.

Finally, issue

.. code-block:: bash

   git push

to push the committed files to Bitbucket.

Git comes with tons of useful commands and being proficient at Git is generally a very useful skill to have. The basics of Git are nicely covered in the first three chapters of the Git book: https://git-scm.com/book/en/v2.

Creating a branch
^^^^^^^^^^^^^^^^^^

When contributing a fix or a new feature, it is best to start with a pull request.

Start by creating a new branch, e.g.

.. code-block:: bash

   git checkout -b demod-fix

Once you have implemented and committed your changes, you can push to Bitbucket by issuing

.. code-block:: bash

   git push --set-upstream origin demod-fix

Next, go to the repository's Bitbucket page (https://bitbucket.org/sron_earth/teds), go to ``Pull requests`` and ``Create pull request``. The default ``Title`` generated from the commit message and an empty ``Description`` field should be fine although both can be edited. Select a reviewer (or leave empty if there is no obvious choice) and tick the ``Delete branch`` field unless you want to keep the branch around and manually delete it later.

Testing out a branch
^^^^^^^^^^^^^^^^^^^^^

In order to test new code written by another developer, first run

.. code-block:: bash

   git pull

to update (merge) all branches. Alternatively, run

.. code-block:: bash

   git fetch

if you want to download the changes but not yet apply them to your local branches. Then issue

.. code-block:: bash

   git checkout new-feature

where ``new-feature`` is the name of the new branch. Next, recompile and rerun the calculation or whichever way you need to test the new feature. If everything works as expected, switch back to the master branch or some other branch you were working on:

.. code-block:: bash

   git checkout master

You may now delete the local branch

.. code-block:: bash

   git branch -d new-feature

in order to clean up but it's not strictly necessary.

If a pull request was created for the ``new-feature`` branch on Bitbucket and you were selected as a reviewer you can now approve the pull request (find the Approve button). The branch is then automatically merged to the master branch and deleted. If it was not set to be automatically deleted upon merge, it may be manually deleted later. Once the feature branch has been merged to the master branch, either automatically by the pull request or manually by the branch creator or anyone with write access, run

.. code-block:: bash

   git pull

again to update the master branch so that it includes the new feature.


Python
-------


Style guide
^^^^^^^^^^^

This project follows the PEP 8 style guide which is universally adopted by most Python projects. It is also the standard used for the Python standard library development and is described in full here: https://peps.python.org/pep-0008/.

In order to see if your code conforms to the standard, configure your editor to highlights parts of the code that do not conform or use an external tool to do so. One such tool is ``flake8`` which is part of ``requirements.txt`` (the executable is in your path if the virtual environment is activated). You can test the correctness of a source file by running ``flake8 file.py``. This tool compares the source file(s) to a set of rules defined by PEP 8 and generates a report per source file. If everything conforms to the standard there should be no output. Writing readable code is important for *i)* reducing the likelihood of future bugs and *ii)* reducing the time it takes for someone (including yourself) to read and contribute to the code.

``flake8_check.sh``, found in the root source directory, checks all Python source files with some exceptions listed in the script. It is run as one of the steps in the Bitbucket pipeline when new code is pushed. If there are style errors in any of the source files the pipeline fails.

In addition to PEP 8, here are some additional rules specific to the TEDS project:

* Do not use non-ASCII symbols.
* Do not use an empty class for the purpose of amending it across functions. Use a dictionary or, better, inherit from an existing class.
* Do not commit commented out code. The code should either be there or not. If not sure, just don't commit that part at all and come back later.
* Do not commit things like TODO lists unless you are convinced they are informative for everybody. You can keep a TODO comment in the code but there is no need to commit it.
* Group imports as follows: standard system libraries followed by local libraries. Separate the groups by blank lines. Within a group, sort the imports alphabetically. This also means that ``from ...`` should come before ``import ...``. Do not import multiple things on one line.


Type hints
^^^^^^^^^^

*coming soon..*


Doc-strings
^^^^^^^^^^^

*coming soon..*


Regression tests
^^^^^^^^^^^^^^^^

When contributing a new feature or fix, it is important not to break anything in other parts of the code. To make sure that previously developed and tested software still performs after a change - in other words that there has not been a regression -  we run *regression tests* before every commit. If any of the tests fail, the conflict must be resolved so that all tests pass.

Besides running tests in your development environment, the test suite is run automatically in a *runner* with each push to Bitbucket. The status of the latest run of the test suite is seen at the repository's overview page: https://bitbucket.org/sron_earth/teds (on the right side you will find something like "Pipeline # for master"). Detailed logs of all tests are found by clicking on Pipelines on the left side of the page.

A new piece of code or a bug fix typically warrants a new test or amendments to existing tests. It is thus normal for tests to keep growing over time and sometimes even exceed the amount of normal code.


C++
---

Style guide
^^^^^^^^^^^

In TEDS, we make use of two C++ linters of which only one is mandatory. Code linting is an automated process that checks code syntax and readability by comparing it to a set of rules. It's basically the same thing as what ``flake8`` does for Python.

The first tool, called ``clang-format``, checks for formatting violations. You can run it by issuing

.. code-block:: bash

   clang-format file.cpp | diff -u file.cpp -

on a source file or

.. code-block:: bash

   bash clang_format_check.py <clang-format>

in the root source directory where ``<clang-format>`` is the ``clang-format`` executable. If the script returns a diff then there are source code formatting errors which should be resolved before committing.

The second tool, called ``clang-tidy``, is a Clang based C++ linter for diagnosing style violations, interface misuse, and violations of best practices. Just like with regression tests, before comitting, it would be good to run ``clang-tidy`` on the source code but it is not enforced at the moment because it will require some effort to make the code fully compliant.

CMake has built-in support for ``clang-tidy`` so all you need to do is run

.. code-block:: bash

   cmake -DCMAKE_CXX_CLANG_TIDY=clang-tidy .

in the build directory and recompile. You can keep this on but if it noticeably slows down compilation you might want to turn it off with

.. code-block:: bash

   cmake -U CMAKE_CXX_CLANG_TIDY .

``clang-format`` is run automatically along with regression tests at each push to the repository whereas ``clang-tidy`` is not. We leave it up to the developer to run ``clang-tidy`` and inspect its output manually for now.


Coding rules
^^^^^^^^^^^^

The general rule is to follow the C++20 standard. Other than that we don't list the rules in detail because ``clang-tidy`` and ``clang-format`` are already quite exhaustive. If those pass then normally the code is correctly formatted.

That said, here is a small selection of rules we want to draw the contributor's attention to:

* The line limit is 80 characters.
* Never use ``use namespace``.
* Use camel case for function and class names and underscores otherwise.
* When writing comments follow the rules of English grammar. Start all comments, if possible, with capitalization. If the comment is one or more whole sentences use normal punctuation. However, if the comment is a single sentence that fits into one line, do not end with a period. Do not end non-sentences with a period.
* Use spaces in argument lists and with most binary operations.
* Always use signed integers - ideally the default ``int`` - over unsigned ones unless there is a compelling reason to do otherwise like if you read or write a NetCDF4 variable that is defined to be unsigned. When using integer types other than the default one, use fixed-width ones. For example, if you need to represent a value larger than 2^31, prefer the 64-bit type ``int64_t``.

Most of those rules are already covered by the aforementioned code linters. For further tips on the best practices of C++ coding, here is an excellent source: https://isocpp.github.io/CppCoreGuidelines/CppCoreGuidelines


Regression tests
^^^^^^^^^^^^^^^^

The C++ tests are written using the Catch2 testing framework. Catch2 needs to be separately installed if not already present on your system. Tests can then be run with

.. code-block:: bash

   cmake --build . --target test

in either the IM or L1A-L1B processor build directory. Each test is an executable in ``tests`` in the build directory so you can also run them manually one by one.


Code coverage
^^^^^^^^^^^^^^^

Code coverage can be inspected by enabling the ``INCLUDE_COVERAGE`` CMake variable. This only works with the GNU compiler and requires LCOV to be installed. The flags defined by ``COVERAGE_FLAGS`` are then appended to compilation flags (the default ``–coverage`` should be fine). If you then recompile and run the tests it shows you the overall coverage rate for lines and functions in the form of a detailed HTML report. Here is an example
of how to programmatically turn on the coverage flag, run the coverage, and
then turn it off again:

.. code-block:: bash

   cmake -DINCLUDE_COVERAGE=ON .
   cmake --build . --target coverage
   cmake -DINCLUDE_COVERAGE=OFF .

From the report you can see the percentage of code lines covered by tests and also a breakdown per source file. Ideally, we should strive for 100% code coverage but in practice that's rarely achieved. Getting to 90% is already pretty good.


Bitbucket pipelines
-------------------

All tests (e.g. code analysis and regression) are run in a Docker container each time new code is pushed to the repository. As a developer, there is normally no need to run the container yourself. You might wish to do so, however, if the regression tests pass on your computer but fail in the runner. Then entering the container allows you to debug the issue in the exact same environment as where the tests are run.

The recipe for how the Docker image is built is found in ``CI/docker_image/Dockerfile``. You can build it yourself, if you wish, by issuing

.. code-block:: bash

   cd CI/docker_image
   sudo docker build -t tango .

or just pull the latest copy of the image from Docker Hub:

.. code-block:: bash

   sudo docker pull raullaasner/tango

The image presents a minimal environment, based on Ubuntu 24.04, with all the TEDS prerequisites installed. You can generate a container from the image and enter it via

.. code-block:: bash

   sudo docker run -it --rm raullaasner/tango

When done, issue ``Ctrl+D`` to exit and delete the container.

Commands that are run inside the container each time new code is pushed are found in ``bitbucket-pipelines.yml`` in the root source directory. Those steps constitute the so-called pipeline. You can see the status of each pipeline at https://bitbucket.org/sron_earth/teds/pipelines. If the pipeline succeeded then on the main repository page, https://bitbucket.org/sron_earth/teds, you can find a green tick mark (usually lower right corner of the page). If the latest pipeline failed then there is a red cross mark. That is a signal to other developers and users that there could be issues with the code and they should not use the most recent version until the issues are resolved.


Contributing to this document
-------------------------------

This documentation is hosted by Read the Docs service and is generated using the Sphinx documentation tool. The markup language used for writing the documentation is called reStructuredText. It is advisable to work through a reStructuredText tutorial and look at the ``rst`` files that make up this documentation before contributing.

When making changes to the documentation, you can view the result by running

.. code-block:: bash

   make html

in the root directory and opening ``build/html/index.html`` in a web browser. When done editing, commit and push to the repository. Read the Docs service will automatically pick up the changes and update https://teds.rtfd.io/ within minutes.


Debugging with GDB
-------------------

A debugger is a tool to run the target program under controlled conditions that allow the programmer to track its operations step by step and monitor changes in computer resources. It can give you more control in pinpointing the source of an unexpected state of the program (e.g. the calculation terminates early or finishes but yields incorrect results) compared to running the program normally (with either release or debug flags). The only requirement for running the GNU Debugger (GDB) is to include the ``-ggdb`` compiler flag. There are many tutorials about the GDB out there so we only list a few example commands here:

- Start the GDB with the IM or L1B executable as an argument and then run with a configuration file as input:

  .. code-block:: bash

     gdb ./tango_l1b.x
     run l1b.yaml

- Insert a breakpoint to monitor line 6 in ``file.cpp``,

  .. code-block:: bash

     break file.cpp:6

  or insert a breakpoint to monitor a function call:

  .. code-block:: bash

     break my_func

- Execute the next program line and step into any function calls in the line,

  .. code-block:: bash

     step

  or step *over* any function calls on the line:

  .. code-block:: bash

     next

- Continue running the program:

  .. code-block:: bash

     continue

- Delete a specified breakpoint:

  .. code-block:: bash

     delete

- Show information about all declared breakpoints:

  .. code-block:: bash

     info breakpoints

- See the value of a variable in the current state of the program:

  .. code-block:: bash

     print ckd

  If ``ckd`` is a class instance the output can be narrowed by specifying a member of that class:

  .. code-block:: bash

     print ckd.n_act

- Display a stack trace of the function calls that lead to a segmentation fault,

  .. code-block:: bash

     backtrace

  or use

  .. code-block:: bash

     where

  which is same as ``backtrace`` but you can use it while you’re still in the middle of the program.

- Run until the current function is finished:

  .. code-block:: bash

     finish

Performance profiling with Perf
--------------------------------

This section is optional and only useful if you're optimizing the C++ code. The IM and L1A-L1B processor are applications of high performance computing with a focus on translating scientific equations into code and optimizing it for speed and memory. While both are important, typically more time is spent on the speed (timing) analysis which means identifying hotspots in the code and attempting to improve the performance in those regions.

The most basic form of timings analysis is looking at the total time it takes for a calculation to run or looking at the timings of individual components as seen in the output of the code. For a more in-depth understanding of where the bottlenecks occur it is better to use a profiling tool. This section describes how to use the Perf tool.

Perf is a performance analyzing tool that ships with the Linux kernel. It can measure different types of events, the most common ones being software events and hardware events. Examples of software events include the CPU clock and page faults while hardware events refer to the number of cycles, instructions retired, L1 cache misses, and many others. It is recommended to work through a Perf tutorial for a full understanding of its capabilities. To get a list of all supported events issue

.. code-block:: bash

   perf list

In this section, we present a few example commands to get you started with using Perf on the L1A-L1B processor. First, recompile the code using normal release flags plus the ``-ggdb`` flag. The ``perf stat`` command keeps a running count of events during execution and presents a summary at the end of the calculation. For instance, running

.. code-block:: bash

   perf stat -e cycles,instructions,cache-references,cache-misses,L1-dcache-loads,\
   L1-dcache-load-misses,branches,branch-misses tango_l1b.x im.yaml

where the ``-e`` flag specifies which events are measured, will output something like::

  30,284,358,788 cycles:u                (62.47%)
  69,653,498,559 instructions:u          (62.49%) # 2.30 insn per cycle
   2,389,054,378 cache-references:u      (62.58%)
     259,681,295 cache-misses:u          (62.55%) # 10.870 % of all cache refs
  12,353,196,649 branches:u              (62.50%)
      64,626,303 branch-misses:u         (62.49%) # 0.52% of all branches
  25,217,043,707 L1-dcache-loads:u       (62.49%)
   1,120,458,843 L1-dcache-load-misses:u (62.46%) # 4.44% of all L1-dcache accesses
     9.604108489 seconds time elapsed
     7.241902000 seconds user
     2.172071000 seconds sys

The absolute number of events such as CPU cycles or instructions are usually not very meaningful. Ratios such as instructions per cycle (IPC) or the number of CPU cache misses vs all cache access attempts are a better measure of performance. A good value for IPC depends on the processor.

In order to identify the hotspots, i.e. to measure events attributed to a specific function or line of code, Perf does event-based sampling which is a statistical process. This means that not every event is explicitly counted. Instead, a sample is registered after a certain number of CPU cycles have passsed. The number of events attributed to a section of the code is thus approximate and for a low number of events care must be taken in interpreting the results. Identification of hotspots is a 2-step process. First collect the samples with ``perf record``:

.. code-block:: bash

   perf record -e cycles,instructions,cache-references,cache-misses,branches,\
   branch-misses,L1-dcache-loads,L1-dcache-load-misses tango_l1b im.yaml

Then analyze the results using ``perf report``:

.. code-block:: bash

   perf report --percent-limit 0.1 --dsos tango_l1b --stdio --fields overhead,sample,symbol

This displays the number of various events at the function level. In order to analyze the events with source code line numbers, issue

.. code-block:: bash

   perf report --percent-limit 0.1 --dsos tango_l1b --stdio --fields overhead,sample,srcline
