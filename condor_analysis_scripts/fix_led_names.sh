#!/usr/bin/env bash

for file in *_result.root ; do mv $file ${file//ON/On} ; done
