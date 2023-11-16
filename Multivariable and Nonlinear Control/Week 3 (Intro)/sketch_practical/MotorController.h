class MotorController_c {
public:
  MotorController_c() {}

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

  // void function does not return anything
  void SetupMotorController(int pwm_pin, int in_1, int in_2) {
    PWM_PIN = pwm_pin; // these assigned variables are whatever is inputted to the func e.g pwm_Pin inptu may be 5
    IN_1_PIN = in_1;
    IN_2_PIN = in_2;
    pinMode(PWM_PIN, OUTPUT); // set the pins as output
    pinMode(IN_1_PIN, OUTPUT);
    pinMode(IN_2_PIN, OUTPUT);
  }

  void SetMotorPower(int dir, int pwmVal) {
    analogWrite(PWM_PIN, pwmVal); // pwm_pin is pin number on arduino. analogWrite helps control the speed of a motor
    // pwm stands for pulse with modulation, the amount of time for which the signal is on
    if (dir == 1) {
      digitalWrite(IN_1_PIN, HIGH); // HIGH is 1 AND LOW is 0
      digitalWrite(IN_2_PIN, LOW);
    } else if (dir == -1) {
      digitalWrite(IN_1_PIN, LOW);
      digitalWrite(IN_2_PIN, HIGH);
    } else {
      digitalWrite(IN_1_PIN, LOW);
      digitalWrite(IN_2_PIN, LOW);
    }
  }

  void ControlLoop(int encoder_count) {
    e = target_counts - encoder_count; 
    if (e < 0) u = -50; // if encoder counts> target counts
    else if (e > 0) u = 50;
    // sets motor power as abs of control signal
    u_amplitude = abs(u);
    // motor direction based on the sign of control signal
    u_sign = 1;
    if (u < 0) {
      u_sign = -1;
    }
    SetMotorPower(u_sign, u_amplitude); // now remember SetMotorPower takes in the dir (u_sign) and pwmVal (motor power or u_amplitude)

    // calculate d_time
    current_T = micros();
    delta_T = (current_T - previous_T) / 1e6;
    previous_T = current_T;
    // print target and position to see the power
    interval_count = interval_count + 1;
    if (interval_count >= print_interval) {
      interval_count = 0;
      average_delta_T = (micros() - interval_start) / (1e6 * print_interval);
      interval_start = micros();
      Serial.print("target_counts:");
      Serial.print(target_counts);
      Serial.print(", ");
      Serial.println();
    }
  }

  void SetTargetCounts(int counts) {
    target_counts = counts; // this func assigns the taken in value (counts) to the vatiable target_counts
  }
};
