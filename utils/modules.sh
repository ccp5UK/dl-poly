#!/usr/bin/env bash
#
#

[[ -d /opt/modules ]] && module use /opt/modules || echo "no aditional modules"
[[ -d /opt/easybuild/modules/all ]] && module use /opt/easybuild/modules/all || echo "no easybuild modules"
