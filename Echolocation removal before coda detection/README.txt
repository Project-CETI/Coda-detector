## 🐋 Sperm Whale Click Elimination and Coda Detection (MATLAB)

This MATLAB program detects and removes **sperm whale echolocation clicks** prior to applying a **coda detector**. By eliminating high-SNR echolocation clicks before detection, the system improves the **trade-off between detection probability and false alarm rate** for codas.

---

### ⚙️ Requirements

- **MATLAB** with the following add-on:
  - **MinGW-w64 C/C++/Fortran Compiler**
- Compilation (only once required):
  Open the MATLAB Command Window and run:
  ```matlab
  mex mcc4mot.c
  ```

---

### ▶️ How to Run

1. Open the file `Click_elimination_and_coda_detection.m` in MATLAB.
2. Press the **Run** button.
3. A dialog box will appear — **select the folder** containing your recorded audio files.
4. An interactive figure will open, displaying the recorded signal.
5. Press **any key** to activate the amplitude threshold selection.
6. **Click on the plot** to choose the **minimum allowed amplitude** for the echolocation click detector.
7. A visualization of the **click train detection** process will follow.
   - All detected echolocation clicks will be **excluded** from the subsequent **coda detection** step to reduce false positives.
8. A final dialog box will ask whether to apply a **constrained solution** for coda detection:
   - **Recommended** when codas are immersed in strong echolocation clicks (i.e., high SNR environments).

---

### 📊 Output

- A figure showing the recorded signal with **detected codas marked in different colors**.
- Detection results are automatically saved to a `.csv` file.
- The **arrival times of individual clicks** within each detected coda are also logged.
