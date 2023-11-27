class MotorController_c {

  public:
  //Constructor, must exist.
    MotorController_c() {

    }
  // variables initialised to zero.
  // Should de set by calling setupMotorcontroller method
    int PWM_PIN = 0;
    int IN_1_PIN = 0;
    int IN_2_PIN = 0;

    int u = 0;
    int u_amplitude = 0;
    int u_sign = 0;
    float target_counts = 0;

    int e = 0;
   
    long current_T = 0;
    long previous_T = 0;
    int print_interval = 100;
    int interval_count = 0;
    float average_delta_T = 0;
    int interval_start = 0;
    float delta_T = 0;

    void SetupMotorController(int pwm_pin, int in_1, int in_2){
      // assign local variables to MotorController_c
      // class attributes
      PWM_PIN = pwm_pin; // for example motor1 passes 3 values for the pins and these are set to 0 
      IN_1_PIN = in_1;
      IN_2_PIN = in_2;
      // set up input and ouput pins
      pinMode (PWM_PIN, OUTPUT);
      pinMode (IN_1_PIN, OUTPUT);
      pinMode (IN_2_PIN, OUTPUT);
    }

    void SetMotorPower (int dir, int pwmVal) {
      //set the pwm pin to the appropriate level
      analogWrite(PWM_PIN, pwmVal);
      //set direction
      if (dir == 1) {
        digitalWrite(IN_1_PIN, HIGH);
        digitalWrite(IN_2_PIN, LOW);
      } else if (dir == -1) {
        digitalWrite(IN_1_PIN, LOW);
        digitalWrite (IN_2_PIN, HIGH);
      } else {
        digitalWrite(IN_1_PIN, LOW);
        digitalWrite (IN_2_PIN, LOW);
      }
    }

    void ControlLoop(int encoder_count) {
    // current position error
      e = target_counts - encoder_count;  // difference between the target position (target_counts) and the actual position indicated by the encoder (encoder_count)
      // u = 2*e;
      if( e < 0)u = -50;
      else if (e > 0)u = 50;
    // set motor power as absolute of control signal
      u_amplitude = abs(u);
    // motor direction based on sign of control signal
      u_sign = 1;
      if (u < 0) {
        u_sign = -1;
      }
    // call a function to drive the motor at a set direction and power
      SetMotorPower(u_sign, u_amplitude);
     
      current_T = micros();
      delta_T = (current_T - previous_T)/ 1e6;
      previous_T = current_T;
      interval_count = interval_count + 1;
      if (interval_count >= print_interval) {
        interval_count = 0;
        average_delta_T = (micros() - interval_start)/(1e6*print_interval);
        interval_start = micros();
        Serial.print("target_counts:");Serial.print(target_counts);Serial.print(", ");
        Serial.print("enc:");Serial.print(encoder_count);Serial.print(", ");
        Serial.println();
      }

    }

    
    
};