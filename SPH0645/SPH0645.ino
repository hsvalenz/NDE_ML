/* SPH0645 MEMS Microphone Test (Adafruit product #3421)
 *
 * Forum thread with connection details and other info:
 * https://forum.pjrc.com/threads/60599?p=238070&viewfull=1#post238070
 */

#include <Audio.h>
const byte encoder0pinA = 10;//A pin -> the interrupt pin 0
const byte encoder0pinB = 8;//B pin -> the digital pin 4
byte encoder0PinALast;
int duration;//the number of the pulses
boolean Direction;//the rotation direction


// GUItool: begin automatically generated code
AudioInputI2S            i2s1;           //xy=180,111
AudioFilterStateVariable filter1;        //xy=325,101
AudioAmplifier           amp1;           //xy=470,93
AudioAnalyzeFFT1024      fft1024_1;      //xy=616,102
AudioConnection          patchCord1(i2s1, 0, filter1, 0);
AudioConnection          patchCord2(filter1, 2, amp1, 0);
AudioConnection          patchCord3(amp1, fft1024_1);
// GUItool: end automatically generated code


void setup() {
  
  AudioMemory(50);
  filter1.frequency(30); // filter out DC & extremely low frequencies
  amp1.gain(8.5);        // amplify signal to useful range
   //EncoderInit();//Initialize the module

}

void loop() {
  Serial.begin(115200);
//    Serial.print("Pulse:");
//  Serial.println(duration);
//  duration = 20;
// delay(100);
  if (fft1024_1.available()) {
    //sampleCount += 1024;  // FFT1024 processes 1024 samples per call

    Serial.print("FFT: ");
    for (int i = 0; i < 26; i++) {
      float n = fft1024_1.read(i);
      if (n >= 0.001) {
        Serial.print(n, 3);
        Serial.print(" ");
      } else {
        Serial.print("  --  ");
      }
    }
    Serial.println();

   
  }
}
// void EncoderInit() 
// {
//  Direction = true;//default -> Forward
//  pinMode(encoder0pinB,INPUT);
//  attachInterrupt(0, wheelSpeed, CHANGE);
// }

// void wheelSpeed()
// {
//  int Lstate = digitalRead(encoder0pinA);
//  if((encoder0PinALast == LOW) && Lstate==HIGH)
//  {
//  int val = digitalRead(encoder0pinB);
//  if(val == LOW && Direction)
//  {
//  Direction = false; //Reverse
//  }
//  else if(val == HIGH && !Direction)
//  {
//  Direction = true; //Forward
//  }
//  }
//  encoder0PinALast = Lstate;

//  if(!Direction) duration++;
//  else duration--;
// }

