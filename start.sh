#!/bin/bash
set -e

# Start Redis server in the background (no daemonize avoids PID-file permission issues)
redis-server --loglevel warning &

# Wait until Redis accepts connections
echo "Waiting for Redis..."
until redis-cli ping 2>/dev/null | grep -q PONG; do
    sleep 0.5
done
echo "Redis ready."

# Start RQ worker from the App folder (tasks.py resolves task_results.db relative to cwd)
cd /app/App
# The bootstrap fixes its own working directory/import path and verifies both
# queued functions before accepting jobs.
python worker.py &

# Start Streamlit in the foreground so it receives container signals
exec streamlit run app.py --server.address=0.0.0.0
