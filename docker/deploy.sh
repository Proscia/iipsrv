#!/bin/bash
#
# Copyright (c) 2009-2018. Authors: see NOTICE file.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#

# Logs
echo "/tmp/iip.out {"          > /etc/logrotate.d/iip
echo "  copytruncate"                   >> /etc/logrotate.d/iip
echo "  daily"                          >> /etc/logrotate.d/iip
echo "  rotate 14"                      >> /etc/logrotate.d/iip
echo "  compress"                       >> /etc/logrotate.d/iip
echo "  missingok"                      >> /etc/logrotate.d/iip
echo "  create 640 root root"           >> /etc/logrotate.d/iip
echo "  su root root"                   >> /etc/logrotate.d/iip
echo "}"                                >> /etc/logrotate.d/iip

sysctl -w net.core.somaxconn=2048

/tmp/start-iip.sh

touch /tmp/crontab
echo "NB_IIP_PROCESS=$NB_IIP_PROCESS" >> /tmp/crontab
echo "*/1 * * * * /bin/bash /tmp/check-status.sh >> /tmp/cron.out" >> /tmp/crontab
crontab /tmp/crontab
rm /tmp/crontab

service rsyslog restart
service cron restart

tail -F /tmp/iip.out
