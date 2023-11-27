volatile int encoder_count_volatile_motor1 = 0;
  void ENC_A_Interrupt_Motor1 () {
  int b_state = digitalRead (ENC_B_PIN_M1);
    if (b_state == 1) {
    encoder_count_volatile_motor1++;
    } else {
    encoder_count_volatile_motor1--;
    }
  }

volatile int encoder_count_volatile_motor2 = 0;
  void ENC_A_Interrupt_Motor2 () {
    //Serial.println(encoder_count_volatile_motor2);
  int b_state = digitalRead(ENC_B_PIN_M2);
      if (b_state == 1) {
        encoder_count_volatile_motor2++;
      } else {
        encoder_count_volatile_motor2--;
      }
  }