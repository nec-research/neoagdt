Neoantigen Digital Twin for neoantigen-based personalized vaccine design
========================================================================

This project contains the source code for designing a personalized neoantigen-based cancer vaccine.
It uses a digital twin simulation.

Installation
--------------

This project is written in `python3` and can be installed with `pip`.

.. code-block::

    pip3 install .[all]

Documentation
-------------

The documentation for this project can be built with `sphinx`. The necessary dependencies are installed by pip
when installing either the `all` or `docs` optional dependencies.

.. code-block::

    cd docs
    make html

The documentation will then be in the `docs/build/html` folder.

Docker
-------------

To manually build a Docker image (with Docker installed on your machine):

.. code-block::

    docker build --tag neoag-digital-twin:<version> .
