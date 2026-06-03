# Hardware Coil Compression

This repository includes the PCB design for a 4x4 mixing matrix that can be used to perform MRI hardware compression at 7T. 

<p align="center">
  <img src="pcb_design.png" width="400">
</p>


# IDE code
The code to run the compression matrix is included in IDE/4x4_matrix_settings. 

The user will have to insert the correct I2C addresses for DACX578 and MCP23X17:

```cpp
// DAC and I/O Port Expander Addresses
#define DAC1_ADDR 0x4C   // First DAC (0-7 channels)
#define DAC2_ADDR 0x4A   // Second DAC (8-14 channels)
#define IOEXP_ADDR 0x20  // I/O Port Expander
```

Furthermore, in setOutputVoltage supply voltage at the input of the Op Amps needs to be defined. The Op Amp gain is set to 3 but this can also be adjusted based on the resistor values one uses.

```cpp
void setOutputVoltage(int row, int col, float voltage) {
  int index = row * ncol + col; // Flatten 2D to 1D index
  const float ampGain = 3;
  float dacOutputVoltage = voltage / ampGain;
  const float Vcc = 4.07;  // Vcc of Op Amp
  uint16_t dacValue = (uint16_t)round((dacOutputVoltage / Vcc) * 4095);
 ...
```

An example of settings is included in setting1.h.

# Bill of Materials

The complete bill of materials for this project can be found in BOM.xlsx, comprehensive of links to vendors and cost. 

