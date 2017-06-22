// Please include  <string> <vector> <sstream>
class stringmanage{
  public :

    stringmanage();

    ~stringmanage();

    template <class T>
    std::string ntos(T number);

    std::string addc(std::string origin, int digit, bool front =true, std::string character ="0");

    std::string rmc(std::string origin, std::string character = " ");

    std::vector<std::string> sn(int start, int end, int digit);

    template <class X>
    std::string addcton(X number, int digit, bool front = true, std::string character = "0");

  private :
  
    std::string my_string;
    
    std::vector<std::string> v_string;
};

stringmanage::stringmanage(){
  my_string ="";
  v_string = {""};
  v_string.clear();
}

stringmanage::~stringmanage(){
}

template <class T>
std::string stringmanage::ntos(T number){
  std::stringstream ss;
  ss << number;
  my_string = ss.str();
  return my_string;
}

std::string stringmanage::addc(std::string origin, int digit, bool front, std::string character){
  my_string = origin;
  for(int i = 0; i < digit; i++){
    if (front){
      my_string = character + my_string;
    }
    else{
      my_string = my_string + character;
    }
  }
  return my_string;
}

std::string stringmanage::rmc(std::string origin, std::string character){
  std::string withoutspace = "";
  my_string = origin;
  for (size_t a = 0; a < my_string.size(); a++){
    if ( std::string(&my_string[a]).find(" ") == std::string::npos){
      withoutspace.push_back(my_string[a]);
    }
  } 
  my_string = withoutspace;
  return my_string;
}

std::vector<std::string> stringmanage::sn(int start, int end, int digit){
  v_string.clear();
  std::string m_string;
  for (int i = start; i < (end+1); i++){
    m_string=ntos(std::abs(i));
    int mydigit = m_string.size();
    if (digit >= mydigit){
      if (i < 0){
        if(digit >= (mydigit+1)){
          m_string = addc(m_string,digit-mydigit-1);
          v_string.push_back("_" + m_string);
        }
      }
      else{
        m_string = addc(m_string,digit-mydigit);
        v_string.push_back(m_string);
      }
    }
  }
  return v_string;
}

template <class X>
std::string stringmanage::addcton(X number, int digit, bool front, std::string character){
  my_string = ntos(number);
  if (digit > my_string.size()){
    my_string = addc(my_string,digit-my_string.size(),front,character);
  }
  return my_string;
}
