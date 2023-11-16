#define ENC_A_PIN_M1 18

#define ENC_B_PIN_M1 19

#define PWM_PIN_M1 12

#define IN1_PIN_M1 11

#define IN2_PIN_M1 10

#include "MotorController.h"

#include "Encoders.h"

MotorController_c motorController1; // I gues this is just declaring variables and a type is not needed?

void setup() {

  // this statement is important for communication between arduino and computer
  Serial.begin(230400);



  // set up motor 1 ==============

  // insanstiate the motor controller object, so its methods(funcs) can be called, for example SetupMotorController
  motorController1 = MotorController_c();

  motorController1.SetupMotorController(PWM_PIN_M1, IN1_PIN_M1, IN2_PIN_M1);



  // sets pint no 18 to input and so on, code in setup runs oncee
  pinMode(ENC_A_PIN_M1, INPUT);

  pinMode(ENC_B_PIN_M1, INPUT);

  // attach interrupt to pin ENC A PIN M1

  // attachInterrupt takes a digitalPinToInterrupt, a function, and a trigger type as its arguments
  attachInterrupt(digitalPinToInterrupt(ENC_A_PIN_M1), ENC_A_Interrupt_Motor1, RISING);

  // ==============
}



void loop() {

  // call the control loop of each motor controller with the mos fecent encoder count

  // takes in the volatile encoder_count variable from Encoders.h which may increase see conditions, this is compared to the target_count?
  // should target count not be defined and called before the Control loop in both scripts??
  motorController1.ControlLoop(encoder_count_volatile_motor1);

  // set the target angle in counts

  // to rotate once every 4 seconds as in arduino 1000 is 1 second
  // This is not a linearly increasing target – it sets a new target every four seconds due to integer division. I dont understand how
  // Notice the jitters? This is because the motor is overshooting its target. We are not scaling motor power through the motor’s rotation. You need to fixx
  // this later in the unit
  motorController1.SetTargetCounts(1440 * (millis() / (1000 * 4)));
}