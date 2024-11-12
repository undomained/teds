Accessing the code
=====================

The Tango end-to-end simulator is hosted at Bitbucket and can be cloned by issuing

.. code-block:: bash

   git clone git@bitbucket.org:sron_earth/teds.git

If you encounter an error with this command, make sure your SSH key has been uploaded to Bitbucket. Log in at https://bitbucket.org, go to :menuselection:`Personal settings --> SSH keys --> Add key`, pick any label, and copy the contents of your public key, typically found at :file:`~/.ssh/id_rsa.pub`. This needs to be done separately for each machine from which you wish to access the repository.

Also, make sure to join the Tango Slack channel! This is the most efficient way to communicate with everybody involved in the project. Ask any of the developers to provide you access to the Slack channel.


Terms of use
-------------

TEDS is licensed under the 3-clause BSD license found in the LICENSE file in the root directory of the project.


Contributors
-------------

In order to see a list of contributors, clone the repository and issue

.. code-block:: bash

   git shortlog -s -n --all --no-merges

This prints a list of contributors along with a number of commits by each contributor.
