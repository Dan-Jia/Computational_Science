// g++  -std=c++11 -pthread -o d task1.cpp
//./d
// die Ergebnisse:
// 0 Hello test1
// 0 Hello test2
//
// 1 Hello test2
// 1 Hello test1
//
// 2 Hello test2
// 2 Hello test1
//
// 3 Hello test2
// 3 Hello test1
//
// 4 Hello test2
// 4 Hello test1

#include <atomic>
#include <iostream>
#include <thread>

// Implementation des Mutexs
// struct
struct mutex_type {
  std::atomic_flag lock = ATOMIC_FLAG_INIT;
};

// Mutex Funktion
void mutex_lock(mutex_type* m) {
  while (std::atomic_flag_test_and_set(&(m->lock)))
    ;  // warten bis mutex frei
}

void mutex_unlock(mutex_type* m) { std::atomic_flag_clear(&(m->lock)); }

// Implementation der Barriere
// struct
struct barriere_type {
  std::atomic<int> arrive;  // Anzahl der Threads,die zur Barriere kommen.
  std::atomic<int> leave;   // Anzahl der Threads,die von Barrieren verlassen.
  int flag;                 // Anzahl der Threads,die die Barriere fordert.
};

// Initialisierung
void barriere_init(barriere_type* b, int num) {
  b->arrive = 0;
  b->leave = num;
  b->flag = num;
}

// Barriere Funktion
void barriere(barriere_type* b) {
  while (b->leave != b->flag)
    ;
  b->arrive++;  // warten bis alle Threads von Barriere verlassen

  if (b->arrive == b->flag) {
    b->leave = 0;
    b->arrive = 0;
  } else {
    while (b->arrive > 0)
      ;
  }
  b->leave++;
}

// Implementation des Semaphors
struct semaphor_type {
  int flag;   // die maximale Kapazität
  int count;  // die belegte Threads
  struct mutex_type m;
};

void semaphor_init(semaphor_type* s, int num) {
  s->flag = num;
  s->count = 0;
}

void semaphor_acquire(semaphor_type* s) {
  while (true) {
    while (s->count >= s->flag)
      ;  // warten bis Resource nicht voll belegt
    mutex_lock(&(s->m));
    if (s->count < s->flag) {
      s->count++;
      mutex_unlock(&(s->m));
      return;
    }
    mutex_unlock(&(s->m));
  }
}

void semaphor_release(semaphor_type* s) {
  mutex_lock(&(s->m));
  if (s->count > 0) {
    s->count--;
  }
  mutex_unlock(&(s->m));
}

// Test-Funktionen
struct mutex_type mut;
struct barriere_type bar;
struct semaphor_type sem;
int number = 0;

void test1() {
  semaphor_acquire(&sem);
  for (int i = 0; i < 5; i++) {
    while (number == 2) {
    }
    mutex_lock(&mut);
    std::cout << i << " "
              << "Hello test1" << std::endl;
    number++;
    mutex_unlock(&mut);
    barriere(&bar);
  }
  semaphor_release(&sem);
}

void test2() {
  semaphor_acquire(&sem);
  for (int i = 0; i < 5; i++) {
    while (number == 2) {
    }
    mutex_lock(&mut);
    std::cout << i << " "
              << "Hello test2" << std::endl;
    number++;
    mutex_unlock(&mut);
    std::this_thread::sleep_for(std::chrono::seconds(1));
    barriere(&bar);
  }
  semaphor_release(&sem);
}

void test3() {
  semaphor_acquire(&sem);
  for (int i = 0; i < 5; i++) {
    while (number < 2) {
    }
    mutex_lock(&mut);
    std::cout << std::endl;
    number = 0;
    mutex_unlock(&mut);
    barriere(&bar);
  }
  semaphor_release(&sem);
}

// main Funktion
int main() {
  while (true) {
    // Initialisierung der Funktion
    barriere_init(&bar, 3);
    semaphor_init(&sem, 3);

    // Threads creat
    std::thread one(test1);
    std::thread two(test2);
    std::thread three(test3);

    // Threads join
    one.join();
    two.join();
    three.join();
  }
  return 0;
}
