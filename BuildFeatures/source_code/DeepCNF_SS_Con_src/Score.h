#pragma once
#include <cmath>

// ----- definitions ------//
typedef double Score;
typedef double LogScore;

const LogScore LogScore_ZERO = (LogScore) -2e20;
const LogScore LogScore_ZERO_REDUCED = LogScore_ZERO/2;


/* Computes log(x) */
inline LogScore Score_LOG (Score x)
{
  return log(x);
}

/* Computes exp(x) */
inline Score LogScore_EXP (LogScore x)
{
  return exp(x);
}

/* Computes log(exp(x) + 1) */
inline LogScore LogScore_LOOKUP (LogScore x)
{
  if (x > (LogScore) 50) return x;
  if (x < (LogScore) -50) return (LogScore) 0;
  return log(exp(x)+1);
}

/* Computes sum of two numbers in log space */
inline LogScore LogScore_ADD (LogScore x, LogScore y)
{
  if (x < y){ LogScore t = x; x = y; y = t; }
  if (y <= LogScore_ZERO_REDUCED) return x;
  return LogScore_LOOKUP(x-y) + y;
}

/* Computes sum of two numbers in log space */
inline void LogScore_PLUS_EQUALS (LogScore &x, LogScore y)
{
  if (x < y){ LogScore t = x; x = y; y = t; }
  if (y > LogScore_ZERO_REDUCED) 
    x = LogScore_LOOKUP(x-y) + y;
}
