
#define ENC_A_PIN_M1 18 // connects PIN M1 of the motor driver to pin 18 of the Arduino Mega controller
#define ENC_B_PIN_M1 19 // encoder pin B which reads A channel produces pulses, and the B channel provides the direction information.
#define PWM_PIN_M1 12
#define IN1_PIN_M1 11
#define IN2_PIN_M1 10

#define ENC_A_PIN_M2 20
#define ENC_B_PIN_M2 21
#define PWM_PIN_M2 5
#define IN3_PIN_M2 7
#define IN4_PIN_M2 6

#include "MotorController.h"
#include "Encoders.h"

MotorController_c motorController1; // instanstiate to motor objects
MotorController_c motorController2;

void setup() {
 Serial.begin (230400);

    // set up motor 1 ==============
    motorController1 = MotorController_c();
    motorController1.SetupMotorController(PWM_PIN_M1, IN1_PIN_M1, IN2_PIN_M1);


    // set up motor 2 ==============
    motorController2 = MotorController_c();
    motorController2.SetupMotorController(PWM_PIN_M2, IN3_PIN_M2, IN4_PIN_M2);
    

    // why are you definining
    pinMode (ENC_A_PIN_M2, INPUT);
    pinMode (ENC_B_PIN_M2, INPUT);

    pinMode (ENC_A_PIN_M1, INPUT);
    pinMode (ENC_B_PIN_M1, INPUT);



    // attach interrupt to pin ENC A PIN M1
    attachInterrupt(digitalPinToInterrupt(ENC_A_PIN_M1), ENC_A_Interrupt_Motor1, RISING);
    attachInterrupt(digitalPinToInterrupt(ENC_A_PIN_M2), ENC_A_Interrupt_Motor2, RISING);
}

// void loop() {
//   // call the control loop of each motor controller with the mos fecent encoder count
//   motorController1.ControlLoop(encoder_count_volatile_motor1);
//   // n the target angle in counts
//   motorController1.SetTargetCounts((131*16) * (millis() / (1000*4)));

//   // motorController2.ControlLoop(encoder_count_volatile_motor2);
//   // motorController2.SetTargetCounts((131*16) * (millis() / (1000*4)));
// }

void loop() {
  // call the control loop of each motor controller with the most recent encoder count
  // motorController1.ControlLoop(encoder_count_volatile_motor1);
  motorController1.ControlLoop(encoder_count_volatile_motor1);
  motorController2.ControlLoop(encoder_count_volatile_motor2);
 
  // Define an array with your desired target positions
  int targetPositions1[] = {}; // -ve is anticlockwise
  int targetPositions2[] = {-10,-90,80,63,340,-400,-600,430,-50}; // -ve is clockwise

  int numPositions1 = sizeof(targetPositions1) / sizeof(targetPositions1[0]);
  int numPositions2 = sizeof(targetPositions2) / sizeof(targetPositions2[0]);
  // Calculate the time period for each position
  int timePeriod = 400;  // Time period for one complete cycle in milliseconds

  // Calculate the current position index based on time
  unsigned long currentTime = millis();
  int positionIndex1 = (currentTime / timePeriod) % numPositions1;
  int positionIndex2 = (currentTime / timePeriod) % numPositions2;
  // Set the target position based on the current position index
  int targetPosition1 = targetPositions1[positionIndex1];
  int targetPosition2 = targetPositions2[positionIndex2];
  // motorController1.SetTargetCounts(targetPosition);
  motorController1.SetTargetCounts(targetPosition1);
  motorController2.SetTargetCounts(targetPosition2);
  //Serial.println(encoder_count_volatile_motor1);
}