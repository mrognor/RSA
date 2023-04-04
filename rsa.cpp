#include <vector>
#include <iostream>
#include <math.h>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <fstream>

#include <boost/random/random_device.hpp>
#include <boost/random.hpp>
#include <boost/multiprecision/cpp_int.hpp>

#define RSA_1024 1024
#define RSA_2048 2048
#define RSA_4096 4096

typedef boost::multiprecision::cpp_int longInt;

struct RSA_Keys
{
    RSA_Keys(const std::string& pubKeyFile, const std::string& privKeyName)
    {
        ReadPrivKey(privKeyName);
        ReadPubKey(pubKeyFile);
    }

    RSA_Keys(longInt e, longInt d, longInt n) { E = e; N = n; D = d; };

    // Return true if successfully open file. False otherwise
    // Private key will saved with file name fileName
    // Public key will saved with file name fileName.pub
    bool SaveKeysToFile(const std::string& fileName)
    {
        std::ofstream privKeyFile(fileName);
        std::ofstream pubKeyFile(fileName + ".pub");

        if (!pubKeyFile.is_open() || !privKeyFile.is_open())
            return false;
        
        privKeyFile << D << std::endl;
        privKeyFile << N << std::endl;
        pubKeyFile << E << std::endl;
        pubKeyFile << N << std::endl;

        pubKeyFile.close();
        privKeyFile.close();
        return true;
    }

    // Return true if successfully open file. False otherwise
    bool ReadPubKey(const std::string& fileName)
    {
        std::ifstream file(fileName);

        if (!file.is_open())
            return false;

        std::string fileLine;

        getline(file, fileLine);
        E.assign(fileLine);

        getline(file, fileLine);
        N.assign(fileLine);
        file.close();

        return true;
    }

    // Return true if successfully open file. False otherwise
    bool ReadPrivKey(const std::string& fileName)
    {
        std::ifstream file(fileName);

        if (!file.is_open())
            return false;

        std::string fileLine;

        getline(file, fileLine);
        D.assign(fileLine);

        getline(file, fileLine);
        N.assign(fileLine);
        file.close();

        return true;
    }

    longInt E;
    longInt D;
    longInt N;
};

// A table of constants from primes up to 1000. It is needed for a quick check for the non-primes of large numbers
const std::vector<int> first_primes = { 
    2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 
    67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 
    181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 
    307, 311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 419, 421, 431, 433, 
    439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499, 503, 509, 521, 523, 541, 547, 557, 563, 569, 571, 577, 
    587, 593, 599, 601, 607, 613, 617, 619, 631, 641, 643, 647, 653, 659, 661, 673, 677, 683, 691, 701, 709, 719, 
    727, 733, 739, 743, 751, 757, 761, 769, 773, 787, 797, 809, 811, 821, 823, 827, 829, 839, 853, 857, 859, 863, 
    877, 881, 883, 887, 907, 911, 919, 929, 937, 941, 947, 953, 967, 971, 977, 983, 991, 997 };

// The function is the inverse of multiplication modulo. (a × b) % c == 1, a = ?
longInt InverseMod(const longInt& b, const longInt& c)
{
    if (b < 1 or c < 2)
        return -1;
 
    longInt u1 = c;
    longInt u2 = 0;
    longInt v1 = b;
    longInt v2 = 1;
 
    while (v1 != 0)
    {
        longInt q = u1 / v1;
        longInt t1 = u1 - q*v1;
        longInt t2 = u2 - q*v2;
        u1 = v1;
        u2 = v2;
        v1 = t1;
        v2 = t2;
    }
 
    if (u1 == 1)
        return (u2 + c) % c;
    else
        return -1;
}

// Algorithm for checking the number for primal
bool MillerRabin (const longInt& n) 
{
    if (n <= 2)
        return false;

	longInt b = 2;
	for (longInt g; (g = boost::multiprecision::gcd(n, b)) != 1; ++b)
		if (n > g)
			return false;

	longInt p = 0, q = n - 1;
	while ((q & 1) == 0)
		++p,  q >>= 1;

	longInt rem = boost::multiprecision::powm(b, q, n);

	if (rem == 1 || rem == n-1)
		return true;

	for (longInt i = 1; i < p; ++i) 
    {
		rem = (rem * (longInt)1 * rem) % n;
		if (rem == n - 1)  
            return true;
	}
	return false;
}

longInt StringToInt(const std::string& str)
{
    longInt n = 0;
    for (auto it = str.rbegin(); it != str.rend(); ++it)
    {
        n += (unsigned char)*it;
        n <<= 8;
    }
    n >>= 8;
    return n;
}

std::string IntToString(longInt n)
{
    std::string str;

    while(n > 0)
    {
        str += (unsigned char)(n % 256);
        n /= 256;
    }

    return str;
}

// Finds the nearest prime number to the given
// varToSet - The number from which the search starts. At the end of the function, the result will be stored in this variable.
// varToFind - A number that changes at each step and is checked for primarly
// isEndedThisNumberThreads - bool variable to stop 2 threads that finding nearest primal to varToSet
// isEndedDifferentNumberThreads - bool variable to check if different primal number was found befor this number
// step - this variable set direction to find primal and cycle step. Usually use +2 and -2 because even numbers always not primal(except 2)
// mtx and cv is service variables
void FindPrimal(longInt& varToSet, longInt& varToFind, bool& isEndedThisNumberThreads, bool& isEndedDifferentNumberThreads, int step, std::mutex& mtx, std::condition_variable& cv)
{
    // Variable to check division to first primal numbers
    bool canBePrimal = true;

    do
    {
        // Check division to first primal numbers
        canBePrimal = true;
        for (const auto& it : first_primes)
        {
            if (varToFind % it == 0) 
            {
                canBePrimal = false;
                break;
            }
        }

        // 1. If it can be primal by not dividing to one of first primals
        // 2. Checking Fermat's little theorem
        // 3. Using MillerRabin algo to end primarly tests
        if (canBePrimal && boost::multiprecision::powm(longInt(2), varToFind - 1, varToFind) == 1 && MillerRabin(varToFind))
        {
            mtx.lock();
            varToSet = varToFind;
            mtx.unlock();
            break;
        }
        varToFind += step;
    } 
    while(varToSet == -1);

    mtx.lock();
    // Check if it is last thread
    if (isEndedThisNumberThreads && isEndedDifferentNumberThreads)
        cv.notify_one();

    if (!isEndedThisNumberThreads)
        isEndedThisNumberThreads = true;
    mtx.unlock();
}

// Key generation use 2 big primal numbers. primalNumbersSize set numbers size in bits 
RSA_Keys KeyGen(uint64_t primalNumbersSize)
{
    longInt minRand = boost::multiprecision::pow(longInt(2), primalNumbersSize);
    longInt maxRand = boost::multiprecision::pow(longInt(2), primalNumbersSize + 1);
    boost::random::random_device gen;
    boost::random::uniform_int_distribution<longInt> ui(minRand, maxRand);

    // Service variables to multithreading
    std::mutex mtx, cvmtx;
    std::condition_variable cv;
    std::unique_lock<std::mutex> lk(cvmtx);

    // First random primal
    longInt p = -1;
    longInt p1 = ui(gen);
        if (p1 % 2 == 0) ++p1;
        longInt p2 = p1;

    // Second random primal
    longInt q = -1;
    longInt q1 = ui(gen);
        if (q1 % 2 == 0) ++q1;
        longInt q2 = q1;

    bool isEndedFirstThreads = false;
    bool isEndedSecondThreads = false;

    // p1th find nearest primal number to p from right sight
    // p2th find nearest primal number to p from left sight
    // q1th find nearest primal number to q from right sight
    // q2th find nearest primal number to q from left sight
    std::thread p1th(FindPrimal, std::ref(p), std::ref(p1), std::ref(isEndedFirstThreads), std::ref(isEndedSecondThreads), +2, std::ref(mtx), std::ref(cv));
    std::thread p2th(FindPrimal, std::ref(p), std::ref(p2), std::ref(isEndedFirstThreads), std::ref(isEndedSecondThreads), -2, std::ref(mtx), std::ref(cv));
    std::thread q1th(FindPrimal, std::ref(q), std::ref(q1), std::ref(isEndedSecondThreads), std::ref(isEndedFirstThreads), +2, std::ref(mtx), std::ref(cv));
    std::thread q2th(FindPrimal, std::ref(q), std::ref(q2), std::ref(isEndedSecondThreads), std::ref(isEndedFirstThreads), -2, std::ref(mtx), std::ref(cv));
    p1th.detach();
    p2th.detach();
    q1th.detach();
    q2th.detach();

    cv.wait(lk);

    longInt n = p * q;
    longInt fi = (p - 1) * (q - 1);
    longInt e;
    for (longInt i = ui(gen); i < fi; i++)
    {
        if (boost::multiprecision::gcd(i, fi) == 1 && boost::multiprecision::powm(longInt(2), i - 1, i))
        {
            e = i;
            break;
        }
    }
    longInt d = InverseMod(e, fi);

    return RSA_Keys(e, d, n);
}

std::string RSA_Encrypt(const std::string& dataToEncrypt, const RSA_Keys& keys)
{
    return IntToString(boost::multiprecision::powm(StringToInt(dataToEncrypt), keys.E, keys.N));
}

std::string RSA_Decrypt(const std::string& dataToDecrypt, const RSA_Keys& keys)
{
    return IntToString(boost::multiprecision::powm(StringToInt(dataToDecrypt), keys.D, keys.N));
}

void KeyGeneratorBenchmark(uint64_t iterations)
{
    std::cout << "Example: https://www.javamex.com/tutorials/cryptography/rsa_key_length.shtml#:~:text=Apple%20M1%20for-,comparison,-%3A" << std::endl;
    std::cout << "Start: " << (float)clock() / CLOCKS_PER_SEC << " sec" << std::endl;

    clock_t Start = clock();
    clock_t delta, start, max = 0, min = 111111111111;
    for (int i = 0; i < iterations; i++)
    {
        start = clock();
        RSA_Keys keys = KeyGen(RSA_1024);
        delta = clock() - start;
        if (delta < min) min = delta;
        if (delta > max) max = delta;
    }
    clock_t End = clock();

    std::cout << "End: " << (float)clock() / CLOCKS_PER_SEC << " sec" << std::endl;
    std::cout << "Min: " << (float)min / CLOCKS_PER_SEC << " sec" << std::endl;
    std::cout << "Max: " << (float)max / CLOCKS_PER_SEC << " sec" << std::endl;
    std::cout << "Sredn: " << (float)((End - Start) / iterations) / CLOCKS_PER_SEC << " sec" << std::endl;
}

int main()
{    
    // To-Do
    // 1) Check key sizes
    
    // g++ -O3 rsa.cpp -lboost_random
    // -lboost_random-mgw12-mt-x64-1_81

    // You can check this realisation speed
    // KeyGeneratorBenchmark(100);

    // You can load keys from files in struct constructor
    // RSA_Keys keys("keys.pub", "keys");
    
    RSA_Keys keys = KeyGen(RSA_1024);
    keys.SaveKeysToFile("keys");

    keys.ReadPrivKey("keys");
    keys.ReadPubKey("keys.pub");

    std::string message = "Some long message to encrypt it using rsa. Some long message to decrypt it using rsa.";
    std::string encryptedMessage = RSA_Encrypt(message, keys);
    std::string decryptedMessage = RSA_Decrypt(encryptedMessage, keys);
    std::cout << decryptedMessage << std::endl;
}
