volatile int encoder_count_volatile_motor1 = 0;

void ENC_A_Interrupt_Motor1() {
  int b_state = digitalRead(ENC_B_PIN_M1);
  if (b_state == 1) {
    encoder_count_volatile_motor1++;
  } else {
    encoder_count_volatile_motor1--;
  }
}
