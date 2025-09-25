from skyfield.api import load, wgs84
from skyfield.data import hipparcos
import numpy as np
from datetime import datetime, timezone

def get_star_position_angle(observer_lat, observer_lon, observation_time, star_name):
    """
    Calculate the position angle of a star relative to North-South line.
    
    Args:
        observer_lat (float): Observer latitude in degrees (+ for North)
        observer_lon (float): Observer longitude in degrees (+ for East)  
        observation_time (datetime): UTC observation time
        star_name (str): Star name ('Spica', 'Algol', or 'Beta Lyrae')
    
    Returns:
        dict: Contains azimuth (position angle), altitude, and star visibility info
    """
    
    # Load planetary data and timescale
    ts = load.timescale()
    planets = load('de421.bsp')
    earth = planets['earth']
    
    # Define observer location
    observer = earth + wgs84.latlon(observer_lat, observer_lon)
    
    # Convert datetime to Skyfield time
    t = ts.from_datetime(observation_time)
    
    # Star coordinates (J2000.0 epoch)
    stars = {
        'Spica': {'ra_hours': 13.420833, 'dec_degrees': -11.161111},
        'Algol': {'ra_hours': 3.136389, 'dec_degrees': 40.955556}, 
        'Beta Lyrae': {'ra_hours': 18.834722, 'dec_degrees': 33.362778}
    }
    
    if star_name not in stars:
        raise ValueError(f"Star '{star_name}' not found. Available: {list(stars.keys())}")
    
    # Create star object
    star_data = stars[star_name]
    star = load('hip_main.dat').star(
        ra_hours=star_data['ra_hours'],
        dec_degrees=star_data['dec_degrees']
    )
    
    # Calculate star position as seen from observer
    astrometric = observer.at(t).observe(star)
    alt, az, distance = astrometric.apparent().altaz()
    
    # Position angle is the azimuth (measured from North through East)
    position_angle = az.degrees
    altitude = alt.degrees
    
    # Determine if star is visible (above horizon)
    is_visible = altitude > 0
    
    return {
        'star_name': star_name,
        'position_angle_degrees': position_angle,
        'altitude_degrees': altitude,
        'is_visible': is_visible,
        'azimuth_direction': get_compass_direction(position_angle),
        'observer_location': f"{observer_lat:.2f}°N, {observer_lon:.2f}°E",
        'observation_time_utc': observation_time.strftime('%Y-%m-%d %H:%M:%S UTC')
    }

def get_compass_direction(azimuth):
    """Convert azimuth to compass direction."""
    directions = ['N', 'NNE', 'NE', 'ENE', 'E', 'ESE', 'SE', 'SSE', 
                 'S', 'SSW', 'SW', 'WSW', 'W', 'WNW', 'NW', 'NNW']
    index = round(azimuth / 22.5) % 16
    return directions[index]

def calculate_multiple_stars(observer_lat, observer_lon, observation_time):
    """Calculate position angles for all three stars."""
    stars = ['Spica', 'Algol', 'Beta Lyrae']
    results = {}
    
    for star in stars:
        try:
            results[star] = get_star_position_angle(
                observer_lat, observer_lon, observation_time, star
            )
        except Exception as e:
            results[star] = {'error': str(e)}
    
    return results

# Example usage
if __name__ == "__main__":
    # Connecticut coordinates (approximate)
    ct_lat, ct_lon = 41.6, -72.7
    
    # Negev Desert coordinates (approximate)  
    negev_lat, negev_lon = 31.0, 35.0
    
    # Example observation time (tonight at 10 PM UTC)
    obs_time = datetime(2025, 9, 17, 22, 0, 0, tzinfo=timezone.utc)
    
    print("=== STAR POSITION ANGLES ===\n")
    
    # Calculate for Connecticut
    print("CONNECTICUT OBSERVATIONS:")
    ct_results = calculate_multiple_stars(ct_lat, ct_lon, obs_time)
    
    for star_name, data in ct_results.items():
        if 'error' in data:
            print(f"{star_name}: Error - {data['error']}")
        else:
            print(f"{star_name}:")
            print(f"  Position Angle: {data['position_angle_degrees']:.1f}° ({data['azimuth_direction']})")
            print(f"  Altitude: {data['altitude_degrees']:.1f}°")
            print(f"  Visible: {'Yes' if data['is_visible'] else 'No'}")
            print()
    
    print("\nNEGEV DESERT OBSERVATIONS:")
    negev_results = calculate_multiple_stars(negev_lat, negev_lon, obs_time)
    
    for star_name, data in negev_results.items():
        if 'error' in data:
            print(f"{star_name}: Error - {data['error']}")
        else:
            print(f"{star_name}:")
            print(f"  Position Angle: {data['position_angle_degrees']:.1f}° ({data['azimuth_direction']})")
            print(f"  Altitude: {data['altitude_degrees']:.1f}°")
            print(f"  Visible: {'Yes' if data['is_visible'] else 'No'}")
            print()

# Additional utility: Track a star throughout the night
def track_star_motion(star_name, observer_lat, observer_lon, start_time, hours=8, interval_hours=1):
    """Track a star's position angle over time."""
    from datetime import timedelta
    
    times = []
    angles = []
    
    for i in range(int(hours/interval_hours) + 1):
        time_point = start_time + timedelta(hours=i*interval_hours)
        result = get_star_position_angle(observer_lat, observer_lon, time_point, star_name)
        
        if result['is_visible']:
            times.append(time_point.strftime('%H:%M'))
            angles.append(result['position_angle_degrees'])
    
    return times, angles