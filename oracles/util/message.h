#ifndef MESSAGE_H
#define MESSAGE_H

#include <vector>

using namespace std;

class Message {
  vector<int> data;
  
public:
  Message();
  Message(vector<int> d);
  Message(int d);
  int read();  
  int &operator[](int i);
};

#endif // MESSAGE_H
