/*
    Owen Stadlwieser and Martin Rudolf
	Assignment 2: Part 2
	CMPUT274, F2019
*/

#include <Arduino.h>
#include <math.h>

const int serverPin = 13;

/*
    Returns true if arduino is server, false if arduino is client
*/
bool isServer() {
    if (digitalRead(serverPin) == HIGH) {
        return true;
    } else {
        return false;
    }
}

/*
    Compute and return (a*b)%m
    Note: m must be less than 2^31
    Arguments:
        a (uint32_t): The first multiplicant
        b (uint32_t): The second multiplicant
        m (uint32_t): The mod value
    Returns:
        result (uint32_t): (a*b)%m
*/
uint32_t multMod(uint32_t a, uint32_t b, uint32_t m) {
    uint32_t result = 0;
    uint32_t dblVal = a%m;
    uint32_t newB = b;

    // This is the result of working through the worksheet.
    // Notice the extreme similarity with powmod.
    while (newB > 0) {
        if (newB & 1) {
            result = (result + dblVal) % m;
        }
        dblVal = (dblVal << 1) % m;
        newB = (newB >> 1);
    }

    return result;
}


/*
    NOTE: This was modified using our multMod function, but is otherwise the
    function powModFast provided in the lectures.

    Compute and return (a to the power of b) mod m.
      Example: powMod(2, 5, 13) should return 6.
*/
uint32_t powMod(uint32_t a, uint32_t b, uint32_t m) {
    uint32_t result = 1 % m;
    uint32_t sqrVal = a % m;  // stores a^{2^i} values, initially 2^{2^0}
    uint32_t newB = b;

    // See the lecture notes for a description of why this works.
    while (newB > 0) {
        if (newB & 1) {  // evalutates to true iff i'th bit of b is 1 in the i'th iteration
            result = multMod(result, sqrVal, m);
        }
        sqrVal = multMod(sqrVal, sqrVal, m);
        newB = (newB >> 1);
    }

    return result;
}



/** Writes an uint32_t to Serial3, starting from the least-significant
 * and finishing with the most significant byte.
 */
void uint32_to_serial3(uint32_t num) {
    Serial3.write((char) (num >> 0));
    Serial3.write((char) (num >> 8));
    Serial3.write((char) (num >> 16));
    Serial3.write((char) (num >> 24));
}


/** Reads an uint32_t from Serial3, starting from the least-significant
 * and finishing with the most significant byte.
 */
uint32_t uint32_from_serial3() {
    uint32_t num = 0;
    num = num | ((uint32_t) Serial3.read()) << 0;
    num = num | ((uint32_t) Serial3.read()) << 8;
    num = num | ((uint32_t) Serial3.read()) << 16;
    num = num | ((uint32_t) Serial3.read()) << 24;
    return num;
}


/*
    Encrypts using RSA encryption.

    Arguments:
        c (char): The character to be encrypted
        e (uint32_t): The partner's public key
        m (uint32_t): The partner's modulus

    Return:
        The encrypted character (uint32_t)
*/
uint32_t encrypt(char c, uint32_t e, uint32_t m) {
    return powMod(c, e, m);
}


/*
    Decrypts using RSA encryption.

    Arguments:
        x (uint32_t): The communicated integer
        d (uint32_t): The Arduino's private key
        n (uint32_t): The Arduino's modulus

    Returns:
        The decrypted character (char)
*/
char decrypt(uint32_t x, uint32_t d, uint32_t n) {
    return (char) powMod(x, d, n);
}

/*
    Core communication loop
    d, n, e, and m are according to the assignment spec
*/
void communication(uint32_t d, uint32_t n, uint32_t e, uint32_t m) {
    // Consume all early content from Serial3 to prevent garbage communication
    while (Serial3.available()) {
        Serial3.read();
    }

    // Enter the communication loop
    while (true) {
        // Check if the other Arduino sent an encrypted message.
        if (Serial3.available() >= 4) {
            // Read in the next character, decrypt it, and display it
            uint32_t read = uint32_from_serial3();
            Serial.print(decrypt(read, d, n));
        }

        // Check if the user entered a character.
        if (Serial.available() >= 1) {
            char byteRead = Serial.read();
            // Read the character that was typed, echo it to the serial monitor,
            // and then encrypt and transmit it.
            if ((int) byteRead == '\r') {
                // If the user pressed enter, we send both '\r' and '\n'
                Serial.print('\r');
                uint32_to_serial3(encrypt('\r', e, m));
                Serial.print('\n');
                uint32_to_serial3(encrypt('\n', e, m));
            } 
            else {
                Serial.print(byteRead);
                uint32_to_serial3(encrypt(byteRead, e, m));
            }
        }
    }
}


/*
    Performs basic Arduino setup tasks.
*/
void setup() {
    init();
    Serial.begin(9600);
    Serial3.begin(9600);
    Serial.println("Welcome to Arduino Chat!");
}

/* 
	Computes the GCD of two given numbers.

	Arguments: uint32_t a, uint32_t b 
	Returns: uint32_t a 
*/
uint32_t gcd_euclid_fast(uint32_t a, uint32_t b) {
  while (b > 0) {
    a %= b;

    // now swap them
    uint32_t tmp = a;
    a = b;
    b = tmp;
  }
  return a; // b is 0
}

/*
	Finds case where gcd is not 1 (i.e. number is prime).

	Arguments: a, some number whose gcd we want to evaluate
	Return: gcd_nonprime
*/
uint32_t gcd(uint32_t a)
{
    for(uint32_t i = 2; i < sqrt(a); i++)
    {
        uint32_t gcd_nonprime = gcd_euclid_fast(a, i);
        if(gcd_nonprime != 1)
        {
            return gcd_nonprime;
        }
    }
    return 1;
}

/* 
	Generates a random k-bit number. 

	Arguments: uint32_t range (desired k)
	Returns: uint32_t delim (random k-bit number)
*/
uint32_t generate_num(uint32_t range){
    uint32_t delim = 0;
    //Reads last bit of range random numbers and creates a range length bit random number
    for (uint32_t i = 0; i < range; ++i) {
        uint32_t val = analogRead(A1);
        uint32_t sigbit = (val&1);
        delay(5);
        delim = delim | (sigbit << i);
    }
    return delim;
}

/* 
    Description: Using the two primes generated, calculates the Arduino modulus, n.

    Arguments: p, q, n, phi (all uin32_t type and pass-by-reference)
*/
void n_calc(uint32_t& p, uint32_t& q, uint32_t& n, uint32_t& phi){
    uint32_t sprime = 10, lprime = 10;
    uint32_t delim = (1 << 14);
    
    // While either number generated is not prime, adds 1 to it
    while(gcd(sprime) != 1 || gcd(lprime) != 1)
    {
        if(gcd(sprime) != 1)
        {
            sprime = generate_num(14);
            sprime = sprime + (delim);
        }
        if(gcd(lprime) != 1)
        {
            lprime = generate_num(15);
            lprime = lprime  + (delim * 2);
        }
    }
    p = lprime;
    q = sprime;
    n = (p*q);
    phi = (p-1)*(q-1);
}

/*
	Generates a public key to be exchanged in handshake protocol. 
	
	Arguments: uint32_t phi
	Returns: uint32_t e (public key)
*/
uint32_t ard_public_key(uint32_t phi){
    uint32_t e = generate_num(14);
    e = e + (1ul <<14);
    // generate random e co prime with phi
    while (gcd_euclid_fast(e ,phi) != 1)
    {
    	// increment e from random starting point reset if too big
    	e = e + 1;
        if(e == ((1ul <<15)-1));
        {
        	e = generate_num(14);
    		e = e + (1ul <<14);
        }
    }
    return e;
}

/*
	Computes the Arduino private key (the modular inverse of e) using Euclid's 
	Extended Algorithm.
	
	Arguments: uint32_t e, uint32_t phi
	Returns: int32_t s[i-1] (Arduino private key)
*/
uint32_t mod_inverse(uint32_t e, uint32_t phi)
{
	uint32_t r[40];
	int32_t s[40];
	int32_t t[40];
	r[0] = e;
	r[1] = phi;
	s[0] = 1;
	s[1] = 0;
	t[0] = 0;
	t[1] = 1;
	int i = 1;
	
	//Euclid's Extended Algorithm
	while (r[i] > 0)
	{
		uint32_t q = r[i-1]/r[i];
		r[i+1] = r[i-1] - (q* r[i]);
		s[i+1] = s[i-1] - (q * s[i]);
		t[i + 1] = t[i-1] - (t[i]*q);
		i = i+1;
	}
	if(s[i-1] < 0)
	{
		s[i-1] = s[i-1] + r[1];
	}
	return s[i-1];
}

/** Waits for a certain number of bytes on Serial3 or timeout
 * @param nbytes: the number of bytes we want
 * @param timeout: timeout period (ms); specifying a negative number 
 *					turns off timeouts (the funciton waits indefinitely 
 *					if timeouts are turned off).
 * @return True if the required number of bytes has arrived. 	
 */
bool wait_on_serial3(uint8_t nbytes, long timeout){
	unsigned long deadline = millis() + timeout;
	while (Serial3.available()<nbytes && (timeout<0 || millis()<deadline))
	{
		delay(1);
	}
	return Serial3.available() >= nbytes;
}

enum ClientState{
	START, CLIWAITINGFORACK, CLIDATAEXCHANGE
};

enum ServerState{
	LISTEN, WAITINGFORKEY, WAITINGFORACK, WAITINGFORKEY2, DATAEXCHANGE
};

/** Runs server side of communication protocol. 
 *	
 *	@param e: public key of the server arduino
 *	@param n: modulus when acting as a client
 *	@param m: modulus when acting as a server
 *	@param otherkey: public key of the client arduino
 */
void runServer(uint32_t e, uint32_t n, uint32_t&m, uint32_t& otherkey){
	int CR;
	bool bytesarrived;
	int message;

	ServerState state = LISTEN;
	while (state != DATAEXCHANGE){
		if (state == LISTEN){
			if (Serial3.available() > 0 ){
				CR = Serial3.read();
				if(CR == 'C'){
					state = WAITINGFORKEY;
				}
			}
		}

		else if(state == WAITINGFORKEY){
			bytesarrived = wait_on_serial3((uint8_t)8, (long)1000);
			if (bytesarrived){
				otherkey = uint32_from_serial3();
				m = uint32_from_serial3();
				Serial3.write('A');
				// send server keys back
				uint32_to_serial3(e);
				uint32_to_serial3(n);
				state = WAITINGFORACK;
			}
			else{
				state = LISTEN;
			}
		}
		else if (state == WAITINGFORACK){
			bytesarrived = wait_on_serial3((uint8_t)1, (long)1000);
			if (bytesarrived){
				message = Serial3.read();
				if(message == 'A'){
					// key exchange complete
					state = DATAEXCHANGE;
				}
				if(message == 'C'){
					state = WAITINGFORKEY2;
				}
			}
			else{
				state = LISTEN;
			}
		}
		else if (state == WAITINGFORKEY2){
			bytesarrived = wait_on_serial3((uint8_t)8, (long)1000);
			if (bytesarrived){
				otherkey = uint32_from_serial3();
				m = uint32_from_serial3();
				state = WAITINGFORACK;
			}
			else{
				state = LISTEN;
			}
		}
	}
};

/** Runs client side of communication protocol. 
 *	
 *	@param e: public key of the client arduino
 *	@param n: modulus when acting as a server
 *	@param m: modulus when acting as a client
 *	@param otherkey: public key of the server arduino
 */
void runClient(uint32_t e, uint32_t n, uint32_t& m, uint32_t& otherkey){	
	int ACK;
	bool bytesarrived;

	ClientState state = START;
	while (state != CLIDATAEXCHANGE){
		if(state== START){
			Serial3.write('C');
			uint32_to_serial3(e);
			uint32_to_serial3(n);
			state =  CLIWAITINGFORACK;
		}
		else if (state == CLIWAITINGFORACK){
			bytesarrived = wait_on_serial3((uint8_t)9, (long)1000);
			if (bytesarrived){
				ACK = Serial3.read();
				if (ACK == 'A'){
					otherkey = uint32_from_serial3();
					m = uint32_from_serial3();
					Serial3.write('A');
					// key exchange complete
					state = CLIDATAEXCHANGE;	
				}
				else{
					state = START;
				}
			}
			else{
				state = START;
			}
		}
	}
};

/*
    The entry point to our program.
*/
int main() {
    setup();
    uint32_t n, p, q, phi, m, otherkey;
    // pass by reference
    n_calc(p, q, n, phi);
   	uint32_t e = ard_public_key(phi);
   	uint32_t d = mod_inverse(e, phi);

 	bool server = isServer();
 	if (server){
 		runServer(e, n, m, otherkey);
 		communication(d, n, otherkey, m);
 	}
 	else{ 
 		runClient(e, n, m, otherkey);
 		communication(d, n, otherkey, m);
 	}

    Serial.flush();
    Serial3.flush();

    return 0;
}