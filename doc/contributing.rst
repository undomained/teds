Contributing
==============

Git
---

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


Regression tests
----------------

When contributing a new feature or fix, it is important not to break anything in other parts of the code. To make sure that previously developed and tested software still performs after a change - in other words that there has not been a regression -  we run *regression tests* before every commit. If any of the tests fail, the conflict must be resolved so that all tests pass.

The tests are written using Catch2 which is a C++ testing framework. Catch2 needs to be separately installed if not already present on your system. Tests can then be run with

.. code-block:: bash

   cmake --build . --target test

Each test is an executable in ``tests`` in your build directory so you can also run them manually one by one.

In addition to the tests you run in your development environment, the test suite is also run automatically in a *runner* with each push to Bitbucket. The status of the latest run of the test suite is seen at the repository's overview page: https://bitbucket.org/sron_earth/teds (on the right side you will find something like "Pipeline # for master"). Detailed logs of all tests are found by clicking on Pipelines on the left side of the page. Even if the tests passed in your own development environment it is still possible for them to fail in the runner. If so, the cause of the failure(s) should be tracked and resolved.

A new piece of code or a bug fix typically warrants a new test or amendments to existing tests. It is thus normal for tests to keep growing over time and even exceed the amount of normal code.


C++
---

Code linters
^^^^^^^^^^^^

Code linting is an automated process that checks code syntax and readability by comparing it to a set of rules. Writing readable and stable code is important for *i)* reducing the likelihood of future bugs and *ii)* reducing the time it takes for someone (including yourself) to read and contribute to the code. In TEDS, we make use of two linters for that purpose.

``clang-tidy`` is a Clang based C++ linter tool for diagnosing style violations, interface misuse, and violations of best practices. Just like with regression tests, before comitting, you should run ``clang-tidy`` on the source code. CMake has built-in support for that so all you need to do is run

.. code-block:: bash

   cmake -DCMAKE_CXX_CLANG_TIDY=clang-tidy .

in the build directory and recompile. You can keep this on but if it noticeably slows down compilation you might want to turn it off with

.. code-block:: bash

   cmake -U CMAKE_CXX_CLANG_TIDY .

The second code linter used is ``clang-format`` which only checks for formatting violations. You can run it by issuing

.. code-block:: bash

   bash clang_format_check.py <clang-format>

in the root source directory where ``<clang-format>`` is the ``clang-format`` executable. If there are no source code formatting errors the script will report so. Otherwise, the error should be resolved before committing. This requirement does not hold for the tests directory.

``clang-format`` is run automatically along with regression tests at each push to the repository whereas ``clang-tidy`` is not. We leave it up to the developer to run ``clang-tidy`` and inspect its output manually for now.

Code coverage
^^^^^^^^^^^^^^^

Code coverage can be inspected by enabling the ``INCLUDE_COVERAGE`` CMake variable. This only works with the GNU compiler and requires LCOV to be installed. The flags defined by ``COVERAGE_FLAGS`` are then appended to compilation flags (the default ``–coverage`` should be fine). If you then recompile and run the tests it shows you the overall coverage rate for lines and functions in the form of a detailed HTML report. Here is an example
of how to programmatically turn on the coverage flag, run the coverage, and
then turn it off again:

.. code-block:: bash

   cmake -DINCLUDE_COVERAGE=ON .
   cmake --build . --target coverage
   cmake -DINCLUDE_COVERAGE=OFF .

Coding rules
^^^^^^^^^^^^

The general rule is to follow the C++20 standard. Other than that we don't list here every small rule because ``clang-tidy`` and ``clang-format`` are already quite exhaustive. If those pass then normally you're good to go.

That said, here is a small selection of rules we want to draw the contributor's attention to:

* The line limit is 80 characters.
* Never use ``use namespace``.
* Use camel case for function and class names and underscores otherwise.
* When writing comments follow the rules of English grammar. Start all comments, if possible, with capitalization. If the comment is one or more whole sentences use normal punctuation. However, if the comment is a single sentence that fits into one line, do not end with a period. Do not end non-sentences with a period.
* Use spaces in argument lists and with most binary operations.
* Always use signed integers - ideally the default ``int`` - over unsigned ones unless there is a compelling reason to do otherwise, for instance if you read or write a NetCDF4 variable that is defined to be unsigned. When using integer types other than the default one, use fixed-width ones. For example, if you need to represent a value larger than 2^31, prefer the 64-bit type ``int64_t``.

Most of those rules are already covered by the aforementioned code linters. For further tips on the best practices of C++ coding, here is an excellent source: https://isocpp.github.io/CppCoreGuidelines/CppCoreGuidelines


Contributing to this document
-------------------------------

This documentation is hosted by Read the Docs service and is generated using the Sphinx documentation tool. The markup language used for writing the documentation is reStructuredText. It is recommended to go through at least one reStructuredText tutorial and look at the ``rst`` files that make up this documentation before contributing.

In order to contribute, clone the public repository where this documentation is hosted:

.. code-block:: bash

   git clone git@bitbucket.org:sron_earth/teds.git

Sphinx is a Python package and in order to ensure compatibility with the correct Python packages, it is important to work in a Python virtual environment. You can create a virtual environment with

.. code-block:: bash

   python -m venv venv

and activate it with

.. code-block:: bash

   source venv/bin/activate

You will notice that the shell's prompt has changed to remind you that you are in a virtual environment. Any packages installed with the Python ``pip`` command are now part of the current project only. Correct versions of packages that are required for this project are listed in ``requirements.txt``. Install them all by issuing

.. code-block:: bash

   pip install --upgrade pip
   pip install -r requirements.txt

When making changes to the documentation you can view the result by running

.. code-block:: bash

   make html

in the root directory and opening ``build/html/index.html`` in a web browser. When done editing, commit and push to the repository,

.. code-block:: bash

   git push

Read the Docs service will automatically pick up the changes and update https://teds.rtfd.io/ within minutes.


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

     print ckd.dim_fov

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

The IM and L1A-L1B processor are applications of high performance computing with a focus on translating scientific equations into code and optimizing it for speed and memory. While both are important, typically more time is spent on the speed (timing) analysis which means identifying hotspots in the code and attempting to improve the performance in those regions.

The most basic form of timings analysis is looking at the total time it takes for a calculation to run or looking at the timings of individual components as seen in the output of the code. For a more in-depth understanding of where the bottlenecks occur it is better to use a profiling tool. This section describes how to use the Perf tool.

Perf is a performance analyzing tool that ships with the Linux kernel. It can measure different types of events, the most common ones being software events and hardware events. Examples of software events include the CPU clock and page faults while hardware events refer to the number of cycles, instructions retired, L1 cache misses, and many others. It is recommended to work through a Perf tutorial for a full understanding of its capabilities. To get a list of all supported events issue

.. code-block:: bash

   perf list

In this section, we present a few example commands to get you started with using Perf on the L1A-L1B processor. First, recompile the code using normal release flags plus the ``-ggdb`` flag. The ``perf stat`` command keeps a running count of events during execution and presents a summary at the end of the calculation. For instance, running

.. code-block:: bash

   perf stat -e cycles,instructions,cache-references,cache-misses,L1-dcache-loads,\
   L1-dcache-load-misses,branches,branch-misses tango_l1b.x cfg.cfg

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
   branch-misses,L1-dcache-loads,L1-dcache-load-misses tango_l1b cfg.cfg

Then analyze the results using ``perf report``:

.. code-block:: bash

   perf report --percent-limit 0.1 --dsos tango_l1b --stdio --fields overhead,sample,symbol

This displays the number of various events at the function level. In order to analyze the events with source code line numbers, issue

.. code-block:: bash

   perf report --percent-limit 0.1 --dsos tango_l1b --stdio --fields overhead,sample,srcline
