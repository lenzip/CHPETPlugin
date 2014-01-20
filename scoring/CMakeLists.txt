CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
PROJECT(ScoringExample)

FIND_PACKAGE(GGS REQUIRED)

ADD_LIBRARY(ScoringExample SHARED 
              src/HitScoringExample.cpp 
              src/SDScoringExample.cpp
              src/UAScoringExample.cpp
            )