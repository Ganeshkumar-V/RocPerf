# RocPerf: An OpenFOAM-based software tool for internal ballistic simulations of solid rocket motors.

![Simulation GIF](docs/media.gif)

RocPerf is a comprehensive, OpenFOAM-based framework for simulating high-speed, gas-particle multiphase flows, capturing the time-resolved evolution of flow fields over the firing duration of rocket motors.

It's a powerful tool for researchers and engineers to study complex internal ballistics, capture propellant regression, and accurately evaluate key performance metrics like thrust, C*, and Isp.
Check out more information in the github pages: https://ganeshkumar-v.github.io/RocPerf/ 

## 🛠️ Requirements

* **Platform:** Linux 

* **Core:** OpenFOAM (developed and tested on **v2112**. This version is recommended, though it may work on other versions.)

* **Compiler:** C++ (GCC or Intel)

* **Build System:** `make` and OpenFOAM toolchain (wmake, Allwmake)

Note: This repository contains a `bashrc` at the project root to set up the environment. The `platforms/` directory contains the library and the applications build on this project

## 🚀 Quick Start — Compilation

1. Ensure OpenFOAM v2112 is installed and source its environment before building this project.

2. Clone the repository in your preferred parent directory (e.g., your projects folder), then enter the project directory:

    ```sh
    git clone https://github.com/Ganeshkumar-V/RocPerf.git
    cd RocPerf
    ```

3.  Source the project `bashrc` to set up the environment:

    ```sh
    source bashrc
    ```

4.  Compile the libraries, applications, and utilities:

    ```sh
    ./Allwmake
    ```

    (Optional: run `./Allwclean` before building to clean the project.)

## 🤝 Contributing

Contributions are welcome. Please open an issue to describe the change, then fork the repository and create a pull request. Follow the existing code style and respect the project's license (GPL-3.0).

## 📄 License

Distributed under the GNU General Public License v3 (GPL-3.0). See the `LICENSE` file for the full text.

This software is an independent project and is not affiliated with, endorsed by, or associated with OpenCFD Ltd or the OpenFOAM® project. OpenFOAM® is a registered trademark of OpenCFD Ltd.

## 🧑‍🔬 Authors and Acknowledgements

* **Lead Developer:** Ganeshkumar V

* **Developed under the supervision of:** Prof. Dilip Srinivas Sundaram

## 📜 Citation

If you use this software in published work, please cite the associated paper:

Venukumar, Ganeshkumar and Sundaram, Dilip Srinivas, "Computational Study of Propulsive Performance of Frozen Nano-Aluminum and Water (ALICE) Mixtures," *Journal of Propulsion and Power*, vol. 41, no. 3, pp. 330–346, 2025. https://doi.org/10.2514/1.B39541

**Suggested BibTeX:**

```bibtex
@article{doi:10.2514/1.B39541,
author = {Venukumar, Ganeshkumar and Sundaram, Dilip Srinivas},
title = {Computational Study of Propulsive Performance of Frozen Nano-Aluminum and Water (ALICE) Mixtures},
journal = {Journal of Propulsion and Power},
volume = {41},
number = {3},
pages = {330-346},
year = {2025},
doi = {10.2514/1.B39541}
}
```

## 📬 Contact

If you encounter problems, please open an issue on the repository.
