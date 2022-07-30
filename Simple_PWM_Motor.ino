
#define motor_pin 3 //define pin 3 for vibration motor 


void setup() {

  //delay(500)
}

void loop() {
  //turn on motor in pin 3, at PWM = 64
  //analogWrite(3, 64);

  //turn on motor in pin 3, at PWM = 128
  //analogWrite(3, 128);

  //turn on motor in pin 3, at PWM = 255
  analogWrite(3, 255);

  
}
