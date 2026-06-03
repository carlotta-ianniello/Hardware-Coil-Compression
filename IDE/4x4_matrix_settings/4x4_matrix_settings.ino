// Include all Adafruit packages
#include <Adafruit_DACX578.h>
#include <Adafruit_MCP23X17.h>

#include "setting1.h"


Adafruit_DACX578 dac1(12);  // Assuming 12-bit resolution
Adafruit_DACX578 dac2(12);  // Assuming 12-bit resolution
Adafruit_MCP23X17 mcp;

// DAC and I/O Port Expander Addresses
#define DAC1_ADDR 0x4C   // First DAC (0-7 channels)
#define DAC2_ADDR 0x4A   // Second DAC (8-14 channels)
#define IOEXP_ADDR 0x20  // I/O Port Expander

const int nrow = 4;
const int ncol = 4;

// Define Pins
const int attPins[] = { 2, 3, 4, 5, 6, 7 };  // 6-bit Attenuator Control Pins
const int LE_PINS[nrow][ncol] = {
  // LE pins for each attenuator
  { 7,  6,  5,  4},     //{A0, A1, A2, A3}    {ch1,1, ch1,2, ch1,3, ch1,4}
  { 3,  2,  1,  0},     //{A4, A5, A6, A7}    {ch2,1, ch2,2, ch2,3, ch2,4}
  {15, 14, 13, 12},     //{B7, B6, B5, B4}    {ch3,1, ch3,2, ch3,3, ch3,4}
  {11, 10,  9,  8}      //{B3, B2, B1, B0}    {ch4,1, ch4,2, ch4,3, ch4,4}
};


void setOutputVoltage(int row, int col, float voltage) {
  int index = row * ncol + col; // Flatten 2D to 1D index
  const float ampGain = 3;
  float dacOutputVoltage = voltage / ampGain;
  const float Vcc = 4.07;  // Vcc of Op Amp
  uint16_t dacValue = (uint16_t)round((dacOutputVoltage / Vcc) * 4095);

  if (index < 8) {
    dac1.writeAndUpdateChannelValue(index, dacValue);
  } else {
    dac2.writeAndUpdateChannelValue(index - 8, dacValue);
  }

  Serial.print("Set (row=");
  Serial.print(row);
  Serial.print(", col=");
  Serial.print(col);
  Serial.print(") → ");
  Serial.print(voltage);
  Serial.println(" V");
}

// Function to set Attenuation values
void setAttenuation(int row, int col, float attenuation_dB) {

  if (attenuation_dB < 0.0) attenuation_dB = 0.0;
  if (attenuation_dB > 31.5) attenuation_dB = 31.5;

  int value = round(attenuation_dB * 2);

  Serial.print("Attenuator [");
  Serial.print(row);
  Serial.print("][");
  Serial.print(col);
  Serial.print("] → ");
  Serial.print(attenuation_dB);
  Serial.println(" dB");

  Serial.print("Binary: ");
  for (int i = 5; i >= 0; i--) {
    int bitValue = (value >> i) & 1;
    digitalWrite(attPins[i], bitValue);
    Serial.print(bitValue);
  }
  Serial.println();

  mcp.begin_I2C(0x20);

  mcp.digitalWrite(LE_PINS[row][col], HIGH);
  Serial.print("LE pin at I/O output: ");
  for (uint8_t i = 0; i < nrow * ncol; i++) {
    int val = mcp.digitalRead(i);
    Serial.print(val);
    Serial.print(" ");
  }
  delay(1);
  mcp.digitalWrite(LE_PINS[row][col], LOW);

  Serial.println();
}

void setup() {

  Serial.begin(9600);
  Wire.setClock(10000); 

  // Initialize  DAC
  if (!dac1.begin(DAC1_ADDR)) {
    Serial.println("Failed to initialize DAC!");
    while (1)
      ;
  }
  if (!dac2.begin(DAC2_ADDR)) {
    Serial.println("Failed to initialize DAC!");
    while (1)
      ;
  }
  // Initialize I/O Expander
  if (!mcp.begin_I2C(IOEXP_ADDR)) {
    Serial.println("Failed to initialize I/O Port Expander!");
    while (1)
      ;
  }

  // Set pins mode
  for (int i = 0; i < 6; i++) {
    pinMode(attPins[i], OUTPUT);
  }

  for (int r = 0; r < nrow; r++) {
    for (int c = 0; c < nrow; c++) {
      mcp.begin_I2C(0x20);
      mcp.pinMode(LE_PINS[r][c], OUTPUT);
      mcp.digitalWrite(LE_PINS[r][c], LOW);
      Serial.println("LE pins set.");
    }
  }

    // Load setting matrixes from file
  Serial.println(B[0][0]);  // Attenuator settings
  Serial.println(D[0][0]);  // Phase shifter settings
}

void loop() {

  for (int r = 0; r < nrow; r++) {
    for (int c = 0; c < nrow; c++) {
      
      setAttenuation(r, c, B[r][c]);
      setOutputVoltage(r, c, D[r][c]);
      
    }
  }

delay(1000);

}
