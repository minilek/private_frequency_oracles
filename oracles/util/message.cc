#include "message.h"

using namespace std;

Message::Message() { }

Message::Message(vector<int> d) {
  data = d;
}

Message::Message(int d) {
  data = vector<int>(1, d);
}

int Message::read() {
  return data[0];
}
  
int &Message::operator[](int i) {
  return data[i];
}
