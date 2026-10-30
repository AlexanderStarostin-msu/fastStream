#!/bin/bash

#git archive --format=zip --output=repo.zip master

git bundle create repo.bundle master
base64 -w 0 repo.bundle > repo.bundle.b64