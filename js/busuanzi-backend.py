#!/usr/bin/env python3
"""
Optional Backend for Standalone Counter
A simple Flask API to sync visitor statistics across devices

Installation:
    pip install flask flask-cors

Usage:
    python3 busuanzi-backend.py

Then in JavaScript:
    const counter = new StandaloneCounter({
        backend: 'http://localhost:5000/api/track'
    });
"""

from flask import Flask, request, jsonify
from flask_cors import CORS
from datetime import datetime
import json
import os

app = Flask(__name__)
CORS(app)  # Enable CORS for cross-origin requests

# Simple file-based storage (use database in production)
DATA_FILE = 'visitor_stats.json'

def load_stats():
    """Load statistics from file"""
    if os.path.exists(DATA_FILE):
        with open(DATA_FILE, 'r') as f:
            return json.load(f)
    return {
        'site_pv': 0,
        'site_uv': 0,
        'visitors': {},
        'pages': {}
    }

def save_stats(stats):
    """Save statistics to file"""
    with open(DATA_FILE, 'w') as f:
        json.dump(stats, f, indent=2)

@app.route('/api/track', methods=['POST'])
def track():
    """Track visitor and return statistics"""
    try:
        data = request.get_json()
        page = data.get('page', '/')
        visitor_id = request.headers.get('X-Visitor-ID')

        # Load current stats
        stats = load_stats()

        # Update site page views
        stats['site_pv'] += 1

        # Track unique visitor
        if visitor_id and visitor_id not in stats['visitors']:
            stats['visitors'][visitor_id] = {
                'first_visit': datetime.now().isoformat(),
                'visit_count': 0
            }
            stats['site_uv'] = len(stats['visitors'])

        if visitor_id:
            stats['visitors'][visitor_id]['visit_count'] += 1
            stats['visitors'][visitor_id]['last_visit'] = datetime.now().isoformat()

        # Track page statistics
        if page not in stats['pages']:
            stats['pages'][page] = {
                'pv': 0,
                'uv': 0,
                'visitors': []
            }

        stats['pages'][page]['pv'] += 1

        if visitor_id and visitor_id not in stats['pages'][page]['visitors']:
            stats['pages'][page]['visitors'].append(visitor_id)
            stats['pages'][page]['uv'] = len(stats['pages'][page]['visitors'])

        # Save updated stats
        save_stats(stats)

        # Return statistics for this page
        return jsonify({
            'site_pv': stats['site_pv'],
            'site_uv': stats['site_uv'],
            'page_pv': stats['pages'][page]['pv']
        })

    except Exception as e:
        return jsonify({'error': str(e)}), 500

@app.route('/api/stats', methods=['GET'])
def get_stats():
    """Get all statistics"""
    stats = load_stats()
    return jsonify(stats)

@app.route('/api/reset', methods=['POST'])
def reset():
    """Reset all statistics (for testing)"""
    if os.path.exists(DATA_FILE):
        os.remove(DATA_FILE)
    return jsonify({'message': 'Statistics reset'})

@app.route('/')
def index():
    """API info page"""
    return '''
    <h1>Busuanzi Backend API</h1>
    <h2>Endpoints:</h2>
    <ul>
        <li><strong>POST /api/track</strong> - Track visitor (body: {page: "/path", timestamp: 123})</li>
        <li><strong>GET /api/stats</strong> - Get all statistics</li>
        <li><strong>POST /api/reset</strong> - Reset statistics</li>
    </ul>
    <h2>Usage:</h2>
    <pre>
const counter = new StandaloneCounter({
    backend: 'http://localhost:5000/api/track'
});
counter.init();
    </pre>
    '''

if __name__ == '__main__':
    print('=' * 50)
    print('Busuanzi Backend API')
    print('=' * 50)
    print(f'Server running at: http://localhost:5000')
    print(f'Track endpoint: http://localhost:5000/api/track')
    print(f'Stats endpoint: http://localhost:5000/api/stats')
    print('=' * 50)
    app.run(debug=True, host='0.0.0.0', port=5000)
