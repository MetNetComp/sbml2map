"""Display-only ordering shared by the SVG and Escher coverage controls."""
import math


def reaction_display_order(records, pathway_priorities=None):
    """Red > orange > other; absolute sampled median and ID resolve ties."""
    priorities = pathway_priorities or {}
    def key(record):
        rid = record['reaction_id']
        priority = max(2 if record.get('is_major_pathway', False) else 0, priorities.get(rid, 0))
        median = float(record.get('median', 0))
        if not math.isfinite(median):
            raise ValueError('Nonfinite sampled median: ' + rid)
        return -priority, -abs(median), rid
    return [record['reaction_id'] for record in sorted(records, key=key)]
